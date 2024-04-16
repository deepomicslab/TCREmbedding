from TCRembedding import get_embedding_instance
import numpy as np

encoder = get_embedding_instance("EmbeddingImRex")

encoder.load_data("data/testdata_ImRex.csv")
encode_result = encoder.embed()

iter_tf_dataset = iter(encode_result)
paired_map_list = []

for item in iter_tf_dataset:
    paired_map, affinity = item
    paired_map_list.append(paired_map.numpy())

encode_result = np.stack(paired_map_list)
print(encode_result.shape)
