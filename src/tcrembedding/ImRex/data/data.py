import pandas as pd

from config import PROJECT_ROOT
from data.control_cdr3_source import ControlCDR3Source
from data.data_source import DataSource
from processing.negative_sampler import add_negatives

DATA_PATH = PROJECT_ROOT / "data/testdata.csv"
allowed_letters = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X']


class DATASource(DataSource):
    """Object holding VDJDB data.
    Contains a pandas DataFrame and dictionary with header names.
    Implements an __iter__ method, and consequently can
    be iterated through via a loop or list comprehension to yield
    the cdr3 and epitope sequences as a tuple, plus a 1.

    Inherits from DataSource object,
    which in turn inherits from Stream object.
    """

    def __init__(
            self,
            filepath=DATA_PATH,
            headers={"cdr3_header": "CDR3", "epitope_header": "Epitope"},
            sep=",",
    ):
        super().__init__()
        self.filepath = filepath
        self.data = pd.read_csv(self.filepath, sep=sep, header=0)
        self.headers = headers

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        for index, row in self.data.iterrows():
            pep1 = row[self.headers["cdr3_header"]]
            pep2 = row[self.headers["epitope_header"]]
            label = row["Affinity"]
            yield (pep1, pep2), label

    def filterseqs(self):
        pattern = '[^' + ''.join(allowed_letters) + ']'
        filtered_df = self.data[~self.data['Epitope'].str.contains(pattern) & ~self.data['CDR3'].str.contains(pattern)]
        self.data = filtered_df.reset_index(drop=True)


    def add_pos_labels(self):
        assert (
                "Affinity" not in self.data.columns
        ), "Dataset already contains class label column Affinity."
        self.data["Affinity"] = 1

    def generate_negatives_from_ref(self, negative_source: ControlCDR3Source):
        """ Generate negative CDR3 epitope sequence pairs by drawing from a negative CDR3 reference set.

        Every epitope in the positive set is matched with a random CDR3 in order to keep the epitope distribution equal between the two classes.
        """
        # sample required number of CDR3 sequences
        amount = self.data.shape[0]

        negative_cdr3_series = (
            negative_source.data[negative_source.headers["cdr3_header"]]
            .sample(n=amount, random_state=42)  # + 3458)
            .reset_index(drop=True)
            .rename(self.headers["cdr3_header"])
        )

        # match with positive epitopes
        negative_df = pd.concat(
            [negative_cdr3_series, self.data[self.headers["epitope_header"]]], axis=1
        )

        # add class labels
        negative_df["Affinity"] = 0

        # merge with positive dataset
        self.data = self.data.append(negative_df).reset_index(drop=True)

        # remove false negatives (and accidental negative duplicates)
        to_do_df = self.data.loc[
            self.data.duplicated(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="last",
            )
        ]
        self.data = self.data.drop_duplicates(
            subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
            keep="first",
        ).reset_index(drop=True)

        # create new negative pairs for any accidental false negatives (and accidental negative duplicates)
        amount = to_do_df.shape[0]
        seed = 42
        while amount > 0:
            seed += 1
            negative_cdr3_series = (
                negative_source.data[negative_source.headers["cdr3_header"]]
                .sample(n=amount, random_state=seed)
                .reset_index(drop=True)
                .rename(self.headers["cdr3_header"])
            )

            # merge with unused epitopes in to_do_df, reset indexing to allow concat
            negative_df = pd.concat(
                [
                    negative_cdr3_series,
                    to_do_df[self.headers["epitope_header"]].reset_index(drop=True),
                ],
                axis=1,
            )

            negative_df["Affinity"] = 0

            self.data = self.data.append(negative_df).reset_index(drop=True)

            to_do_df = self.data.loc[
                self.data.duplicated(
                    subset=[
                        self.headers["cdr3_header"],
                        self.headers["epitope_header"],
                    ],
                    keep="last",
                )
            ]
            amount = to_do_df.shape[0]

            self.data = self.data.drop_duplicates(
                subset=[self.headers["cdr3_header"], self.headers["epitope_header"]],
                keep="first",
            ).reset_index(drop=True)

    def generate_negatives_via_shuffling(
            self, full_dataset_path: str, epitope_ratio: bool = False
    ):
        """Generate negative CDR3-epitope pairs through shuffling and add them to the underlying DataFrame.

        Parameters
        ----------
        full_dataset_path : str
            Path to the entire cdr3-epitope dataset, before splitting into folds, restricting length or downsampling. Used to avoid generating false negatives during shuffling. Should only contain positive values. Will be merged with current train/val dataframe.
            Length trimming = OK
            CV folds =  not OK, in the grouped-kfold setting it does not matter, because when a certain CDR3 is paired with two different epitopes, and they end up in different folds, it's impossible for the CDR3 to be accidentally matched up to the other epitope again, because it's not available for selection. In the normal CV setting it could matter though.
            Downsampling = not OK, a CDR3 could lose representative samples of it being paired with specific epitopes, and could end up being paired with them again as false negatives during shuffling.
            MHC = OK, a few CDR3s occur for both classes, but none of the epitopes do. Consequently it's impossible for a CDR3 to be paired with an epitope that could be a false negative in the full dataset.
            TRAB = OK, none of the CDR3s are identical between TRA and TRB genes. Consequently it's impossible for a CDR3 to be paired with an epitope that could be a false negative in the full dataset.
        epitope_ratio : boolean
            When false, samples an epitope for each CDR3 sequence in the
            proportionally to its occurrence in the other epitope pairs. Does not
            preserve the ratio of positives and negatives within each epitope,
            but does result in every CDR3 sequence having exactly 1 positive and negative.
            When true, samples a set of CDR3 sequences with from the unique list of CDR3s
            for each epitope observation (per epitope), i.e. preserves exact ratio of positives and
            negatives for each epitope, at the expense of some CDR3s appearing more than once
            among the negatives and others only in positives pairs.
        """

        self.data = add_negatives(
            df=self.data,
            full_dataset_path=full_dataset_path,
            epitope_ratio=epitope_ratio,
        )

    def length_filter(
            self,
            min_length_cdr3: int = 10,
            max_length_cdr3: int = 20,
            min_length_epitope: int = 8,
            max_length_epitope: int = 13,
    ):
        self.data = self.data.loc[
            (self.data[self.headers["cdr3_header"]].str.len() >= min_length_cdr3)
            & (self.data[self.headers["cdr3_header"]].str.len() <= max_length_cdr3)
            & (
                    self.data[self.headers["epitope_header"]].str.len()
                    >= min_length_epitope
            )
            & (
                    self.data[self.headers["epitope_header"]].str.len()
                    <= max_length_epitope
            )
            ]
