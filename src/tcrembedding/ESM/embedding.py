import torch
from .esm import pretrained, FastaBatchedDataset, Alphabet

class EmbeddingESM:
    def __init__(self, toks_per_batch=4096, repr_layers=[-1], include=["mean"], truncation_seq_length=1022, nogpu=False):
        self.toks_per_batch = toks_per_batch
        self.repr_layers = repr_layers
        self.include = include
        self.truncation_seq_length = truncation_seq_length
        self.nogpu = nogpu
        self.model = None
        self.alphabet = None

    def load_model(self, model_location):
        self.model, self.alphabet = pretrained.load_model_and_alphabet(model_location)
        self.model.eval()
        if torch.cuda.is_available() and not self.nogpu:
            self.model = self.model.cuda()
            print("Transferred model to GPU")

    def get_fasta_ids(self, fasta_file):
        ids = []
        with open(fasta_file) as f:
            for line in f:
                if line.startswith('>'):
                    id = line.strip().split('>')[1]
                    ids.append(id)
        return ids

    def run(self, model_location, fasta_file):
        self.load_model(model_location)

        smp_ids = self.get_fasta_ids(fasta_file)
        dataset = FastaBatchedDataset.from_file(fasta_file)
        batches = dataset.get_batch_indices(self.toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset, collate_fn=self.alphabet.get_batch_converter(self.truncation_seq_length), batch_sampler=batches
        )

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

            print("Data saved. First embedding:", all_mean_reprs[0])

            return all_mean_reprs
        
if __name__ == "__main__":
    model_location = "esm1b_t33_650M_UR50S"
    fasta_file = "data/IEDB_uniqueTCR_top10_filter.fasta"
    toks_per_batch = 2048
    repr_layers = [33]
    include = ["mean"]
    truncation_seq_length = 1022
    nogpu = True

    extractor = EmbeddingESM(
        toks_per_batch=toks_per_batch,
        repr_layers=repr_layers,
        include=include,
        truncation_seq_length=truncation_seq_length,
        nogpu=nogpu
    )

    extractor.run()