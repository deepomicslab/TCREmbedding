from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingESM")

model_location = "esm1b_t33_650M_UR50S"
fasta_file = "/media/lihe/TCR/project/src/TCRembedding/ESM/data/IEDB_uniqueTCR_top10_filter.fasta"
encoder.toks_per_batch = 2048
encoder.repr_layers = [33]
encoder.include = ["mean"]
encoder.truncation_seq_length = 1022
encoder.nogpu = True

encode_result = encoder.run(model_location, fasta_file)
print(encode_result.shape)