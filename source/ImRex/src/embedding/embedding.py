import numpy as np

from src.bio.feature_builder import CombinedPeptideFeatureBuilder
from src.bio.peptide_feature import parse_features, parse_operator
from src.config import PROJECT_ROOT
from src.data.control_cdr3_source import ControlCDR3Source
from src.data.data import DATASource
from src.processing.data_stream import DataStream
from src.processing.padded_dataset_generator import (
    # augment_pairs,
    padded_dataset_generator,
)

features_list = parse_features("hydrophob,isoelectric,mass,hydrophil,charge")
operator = parse_operator("absdiff")
feature_builder = CombinedPeptideFeatureBuilder(features_list, operator)

class EmbeddingImRex:

    def __init__(self, cdr3_range=(10, 20), create_neg_dataset=False):
        """
        cdr3_range : Tuple[int, int]
        The minimum and maximum desired cdr3 sequence length.
        create_neg_dataset : Union[bool, str]
        Whether to create negatives by shuffling/sampling, by default True.
        NOTE: Should always be set to False when evaluating a dataset that already contains negatives.
        """
        self.cdr3_range = cdr3_range
        self.create_neg_dataset = create_neg_dataset

    def read_csv(self, filepath=r"test_data.csv", negative_data_filepath="negative_data.csv", full_dataset_path="positive_data.csv", cdr3_header="CDR3", epitope_header="Epitope"):
        """
        full_dataset_path : Path
        The entire cdr3-epitope dataset, before splitting into folds, restricting length or downsampling. Used to avoid
        generating false negatives. Should only contain positive values.
        """
        self.filepath = filepath
        self.negative_data_filepath = negative_data_filepath
        self.full_dataset_path = full_dataset_path
        self.cdr3_header = cdr3_header
        self.epitope_header = epitope_header

    def encode(self):
        if self.create_neg_dataset == False:
            data_source = DATASource(
                filepath=self.filepath,
                headers={"cdr3_header": self.cdr3_header, "epitope_header": self.epitope_header},
            )
            data_source.filterseqs()

            data_stream = DataStream(data_source)

            tf_dataset = padded_dataset_generator(
                data_stream=data_stream,
                feature_builder=feature_builder,
                cdr3_range=(self.cdr3_range[0], self.cdr3_range[1]),
                epitope_range=(8, 11),
                neg_shuffle=False
            )

            return tf_dataset
        # repeat without negative reference set
        else:
            data_source = DATASource(
                filepath=self.filepath,
                headers={"cdr3_header": self.cdr3_header, "epitope_header": self.epitope_header},
            )
            data_source.filterseqs()

            if self.create_neg_dataset == "shuffling":
                data_stream = DataStream(data_source)

                tf_dataset = padded_dataset_generator(
                    data_stream=data_stream,
                    feature_builder=feature_builder,
                    cdr3_range=(self.cdr3_range[0], self.cdr3_range[1]),
                    epitope_range=(8, 11),
                    neg_shuffle=True,
                    full_dataset_path=self.full_dataset_path,
                )

                return tf_dataset
            elif self.create_neg_dataset == "sampling":

                negative_source = ControlCDR3Source(
                    filepath=self.negative_data_filepath,
                    min_length=self.cdr3_range[0],
                    max_length=self.cdr3_range[1],
                )

                data_source.generate_negatives_from_ref(negative_source)
                data_stream = DataStream(data_source)

                tf_dataset = padded_dataset_generator(
                    data_stream=data_stream,
                    feature_builder=feature_builder,
                    cdr3_range=(self.cdr3_range[0], self.cdr3_range[1]),
                    epitope_range=(8, 11),
                    neg_shuffle=False
                )

                return tf_dataset

if __name__ == "__main__":

    encoder = EmbeddingImRex()
    encoder.read_csv(filepath=r"D:\TCR\combined_dataset_filtered_length_pos.csv")
    encode_result = encoder.encode()

    iter_tf_dataset = iter(encode_result)
    paired_map_list = []

    for item in iter_tf_dataset:
        paired_map, affinity = item
        paired_map_list.append(paired_map.numpy())

    np.save(f"D:\TCR\combined_dataset\ImRex.npy", paired_map_list)