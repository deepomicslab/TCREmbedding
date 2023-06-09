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

    def __init__(self):
        self.filepath = "testdata.csv"
        self.cdr3_header = "CDR3"
        self.epitope_header = "Epitope"
        self.negative_data_filepath = "negative_data.csv"
        self.cdr3_range = (10, 20)
        self.full_dataset_path = "positive_data.csv"

    def embed(self, filepath, full_dataset_path, create_neg_dataset):
        self.filepath = filepath

        if create_neg_dataset == False:
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

            if create_neg_dataset == "shuffling":
                data_stream = DataStream(data_source)
                self.full_dataset_path = full_dataset_path

                tf_dataset = padded_dataset_generator(
                    data_stream=data_stream,
                    feature_builder=feature_builder,
                    cdr3_range=(self.cdr3_range[0], self.cdr3_range[1]),
                    epitope_range=(8, 11),
                    neg_shuffle=True,
                    full_dataset_path=self.full_dataset_path,
                )

                return tf_dataset
            elif create_neg_dataset == "sampling":

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
