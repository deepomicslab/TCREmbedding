### ESM
# install env: https://github.com/facebookresearch/esm 
# 读入文件是fasta格式，如果是csv文件，可以参考以下命令转换(第一行命令把tsv文件提取第一列，排除表头，存成一个每行是一个tcr序列的文本文件。第二行命令把tcr的文本文件抓换成fasta格式)：
tail -n +2 filename.tsv | awk '{print $1}' > filename.tcr
index=1;for i in `cat filename.tcr` ;do echo '>'$index && echo $i && let index++ ;done  > filename.fasta

# 在ESM目录下执行
python extract_33layer_clustering.py esm1b_t33_650M_UR50S test.fasta ./output_file output.npy --repr_layers 33 --include mean --toks_per_batch 4096

这里esm1b_t33_650M_UR50S是选择的esm模型编号，保持这个选择就可以。
./output_file是输出目录，改成自己的目录
output.npy是输出的numpy文件文件名，改成需要的文件名
--repr_layers 33 --include mean 不用管，是esm设置的参数，这里的意思是要它们transformer模型第33层的输出，而且取平均值当embedding，是github里推荐的做法
--toks_per_batch 4096是batch size。意思是一个batch放多少数据。如果内存够大就可以设置再大些（会跑得快），内存不够就比现在更小。