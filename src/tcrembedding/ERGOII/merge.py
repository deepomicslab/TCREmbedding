import numpy as np

# 加载第一个.npy文件
array1 = np.load('ERGO-II_tcr.npy')

# 加载第二个.npy文件
array2 = np.load('ERGO-II_tcr_1.npy')

# 合并两个数组
merged_array = np.concatenate((array1, array2), axis=0)
print(merged_array.shape)
# 保存合并后的数组为.npy文件
np.save('ERGO-II_tcr.npy', merged_array)
