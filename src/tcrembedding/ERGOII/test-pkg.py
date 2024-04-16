from TCRembedding import get_embedding_instance
from TCRembedding.ERGOII.data_loader import Data_Loader

encoder = get_embedding_instance("EmbeddingERGO")
encoder.tcr_encoding_model = "AE"
encoder.load_model(model_path="Models/version_10ve/checkpoints", args_path="Models/10ve/meta_tags.csv")
data_loader = Data_Loader()
tcrb_list, peptide = data_loader.collate("data/testdata_ERGO-II.csv", "AE")
tcrb_batch, pep_batch = encoder.forward(tcrb_list, peptide)
tcrb_encoding, pep_encoding = encoder.embed(tcrb_batch, pep_batch)

tcr_encode_result = tcrb_encoding.detach().numpy()
pep_encode_result = pep_encoding.detach().numpy()
print(tcrb_encoding.shape)
print(pep_encoding.shape)