import pickle
import numpy as np
import argparse
import pathlib
import os
import torch
from esm import pretrained, FastaBatchedDataset, Alphabet

class EmbeddingESM:
    def __init__(self, model_location, fasta_file, output_dir, output_name="embedding.npy", toks_per_batch=4096, repr_layers=[-1], include=["mean"], truncation_seq_length=1022, nogpu=False):
        self.model_location = model_location
        self.fasta_file = fasta_file
        self.output_dir = pathlib.Path(output_dir)
        self.output_name = output_name
        self.toks_per_batch = toks_per_batch
        self.repr_layers = repr_layers
        self.include = include
        self.truncation_seq_length = truncation_seq_length
        self.nogpu = nogpu
        self.model = None
        self.alphabet = None
        self.output_file = os.path.join(self.output_dir, self.output_name)

    def load_model(self):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(self.model_location)
        self.model.eval()
        if torch.cuda.is_available() and not self.nogpu:
            self.model = self.model.cuda()
            print("Transferred model to GPU")

    def get_fasta_ids(self):
        ids = []
        with open(self.fasta_file) as f:
            for line in f:
                if line.startswith('>'):
                    id = line.strip().split('>')[1]
                    ids.append(id)
        return ids

    def run(self):
        self.load_model()

        smp_ids = self.get_fasta_ids()
        dataset = FastaBatchedDataset.from_file(self.fasta_file)
        batches = dataset.get_batch_indices(self.toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(self.truncation_seq_length), batch_sampler=batches
        )

        self.output_dir.mkdir(parents=True, exist_ok=True)
        repr_layers = [(i + self.model.num_layers + 1) % (self.model.num_layers + 1) for i in self.repr_layers]
        all_embs = {}

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                if torch.cuda.is_available() and not self.nogpu:
                    toks = toks.to(device="cuda", non_blocking=True)

                out = self.model(toks, repr_layers=repr_layers)
                representations = {layer: t.to(device="cpu") for layer, t in out["representations"].items()}
                repre_tensor = representations[repr_layers[0]]

                for i, label in enumerate(labels):
                    truncate_len = min(self.truncation_seq_length, len(strs[i]))
                    if "mean" in self.include:
                        all_embs[label] = repre_tensor[i, 1 : truncate_len + 1].mean(0).clone()

            embs = [all_embs[id] for id in smp_ids]

            all_mean_reprs = torch.stack(embs).cpu().numpy()
            # pickle.dump(all_mean_reprs, file=open(f"{self.output_file}.pkl", 'wb+'))
            # np.save(self.output_file, all_mean_reprs)

            print("Data saved. First embedding:", all_mean_reprs[0])

            return all_mean_reprs
        
if __name__ == "__main__":
    model_location = "esm1b_t33_650M_UR50S"
    fasta_file = "data/IEDB_uniqueTCR_top10_filter.fasta"
    output_dir = "."
    output_name = "example_embeddings.npy"
    toks_per_batch = 2048
    repr_layers = [33]
    include = ["mean"]
    truncation_seq_length = 1022
    nogpu = True

    extractor = EmbeddingESM(
        model_location=model_location,
        fasta_file=fasta_file,
        output_dir=output_dir,
        output_name=output_name,
        toks_per_batch=toks_per_batch,
        repr_layers=repr_layers,
        include=include,
        truncation_seq_length=truncation_seq_length,
        nogpu=nogpu
    )

    extractor.run()