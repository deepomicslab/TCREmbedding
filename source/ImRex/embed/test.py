from embedding import EmbeddingImRex

encoder = EmbeddingImRex()
encode_result = encoder.embed("testdata.csv", "positive_data.csv", "shuffling")
print(type(encode_result))
iter_tf_dataset = iter(encode_result)
for item in iter_tf_dataset:
    tensor1, tensor2 = item
    print(tensor1)
    print(tensor2)
    print("-----------------")
