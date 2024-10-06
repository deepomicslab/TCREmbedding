from TCRembedding import get_embedding_instance

encoder = get_embedding_instance("EmbeddingDeepTCR")
encoder.load_model(train_data_directory='Data/data', model_folder_name="Test_Model", Load_Prev_Data=True)
encoder_result = encoder.embed()
print(encoder_result.shape)
