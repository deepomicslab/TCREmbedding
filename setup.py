import setuptools

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name="tcr-embedding",
    version="1.0.0",
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    package_data={
        # 指定每个嵌入方法的 data 文件夹下的所有文件都应该被包含
        'TCRembedding': ['*/data/*', '*/data/blosum/*', '*/models/*', '*/profiles/*', '*/Data/*', 
                         '*/encode/*', '*/Test_Model/*', '*/Models/*', '*/embedding/*', '*/library/*', 
                         '*/TITAN_model/*',],
    },
)