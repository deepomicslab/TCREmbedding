from h5py import Dataset, Group, File

with File(r'D:\TCR\program\encode\pMTnet\library\h5_file\HLA_antigen_encoder_60.h5', "r") as f:

    model_weight_group = f["model_weights"]
    for key in model_weight_group.keys():
        print(model_weight_group[key], model_weight_group[key].name)

    optimizer_weight_group = f["optimizer_weights"]
    for key in optimizer_weight_group.keys():
        print(optimizer_weight_group[key], optimizer_weight_group[key].name)

with File(r'D:\TCR\program\encode\pMTnet\library\h5_file\HLA_antigen_encoder_60.h5', 'r') as f:
    for k in f.keys():
        if isinstance(f[k], Dataset):
            print(f[k].value)
        else:
            print(f[k].name)