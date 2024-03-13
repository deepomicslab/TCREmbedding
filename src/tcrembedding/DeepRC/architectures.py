import numpy as np
import torch
import torch.nn as nn

class SequenceEmbeddingCNN(nn.Module):
    def __init__(self, n_input_features: int, kernel_size: int = 9, n_kernels: int = 32, n_layers: int = 1):
        """Sequence embedding using 1D-CNN (`h()` in paper)

        See `deeprc/examples/` for examples.

        Parameters
        ----------
        n_input_features : int
            Number of input features per sequence position
        kernel_size : int
            Size of 1D-CNN kernels
        n_kernels : int
            Number of 1D-CNN kernels in each layer
        n_layers : int
            Number of 1D-CNN layers
        """
        super(SequenceEmbeddingCNN, self).__init__()
        self.kernel_size = kernel_size
        self.n_kernels = n_kernels
        self.n_layers = n_layers

        if self.n_layers <= 0:
            raise ValueError(f"Number of layers n_layers must be > 0 but is {self.n_layers}")

        # Set random seed here
        torch.manual_seed(1)

        # CNN layers
        network = []
        for i in range(self.n_layers):
            conv = nn.Conv1d(in_channels=n_input_features, out_channels=self.n_kernels, kernel_size=self.kernel_size,
                             bias=True)
            conv.weight.data.normal_(0.0, np.sqrt(1 / np.prod(conv.weight.shape)))
            network.append(conv)
            network.append(nn.SELU(inplace=True))
            n_input_features = self.n_kernels

        self.network = torch.nn.Sequential(*network)

    def forward(self, inputs, *args, **kwargs):
        """Apply sequence embedding CNN to inputs in NLC format.

        Parameters
        ----------
        inputs: torch.Tensor
            Torch tensor of shape (n_sequences, n_sequence_positions, n_input_features).

        Returns
        ---------
        max_conv_acts: torch.Tensor
            Sequences embedded to tensor of shape (n_sequences, n_kernels)
        """
        inputs = torch.transpose(inputs, 1, 2)  # NLC -> NCL
        # Apply CNN
        conv_acts = self.network(inputs)
        # Take maximum over sequence positions (-> 1 output per kernel per sequence)
        max_conv_acts, _ = conv_acts.max(dim=-1)
        return max_conv_acts


class SequenceEmbeddingLSTM(nn.Module):
    def __init__(self, n_input_features: int, n_lstm_blocks: int = 32, n_layers: int = 1, lstm_kwargs: dict = None):
        """Sequence embedding using LSTM network (`h()` in paper) with `torch.nn.LSTM`

        See `deeprc/examples/` for examples.

        Parameters
        ----------
        n_input_features : int
            Number of input features
        n_lstm_blocks : int
            Number of LSTM blocks in each LSTM layer
        n_layers : int
            Number of LSTM layers
        lstm_kwargs : dict
            Parameters to be passed to `torch.nn.LSTM`
        """
        super(SequenceEmbeddingLSTM, self).__init__()
        self.n_lstm_blocks = n_lstm_blocks
        self.n_layers = n_layers
        if lstm_kwargs is None:
            lstm_kwargs = {}
        self.lstm_kwargs = lstm_kwargs

        if self.n_layers <= 0:
            raise ValueError(f"Number of layers n_layers must be > 0 but is {self.n_layers}")

        # LSTM layers
        network = []
        for i in range(self.n_layers):
            lstm = nn.LSTM(input_size=n_input_features, hidden_size=self.n_lstm_blocks, **lstm_kwargs)
            network.append(lstm)
            n_input_features = self.n_lstm_blocks

        self.network = torch.nn.Sequential(*network)

    def forward(self, inputs, sequence_lengths, *args, **kwargs):
        """Apply sequence embedding LSTM network to inputs in NLC format.

        Parameters
        ----------
        inputs: torch.Tensor
            Torch tensor of shape (n_sequences, n_sequence_positions, n_input_features).

        Returns
        ---------
        max_conv_acts: torch.Tensor
            Sequences embedded to tensor of shape (n_sequences, n_kernels)
        """
        inputs = torch.transpose(inputs, 0, 1)  # NLC -> LNC
        output, (hn, cn) = self.network(inputs)
        output = output[sequence_lengths - 1, torch.arange(output.shape[1], dtype=torch.long)]
        return output


class AttentionNetwork(nn.Module):
    def __init__(self, n_input_features: int, n_layers: int = 2, n_units: int = 32):
        """Attention network (`f()` in paper) as fully connected network.
         Currently only implemented for 1 attention head and query.

        See `deeprc/examples/` for examples.

        Parameters
        ----------
        n_input_features : int
            Number of input features
        n_layers : int
            Number of attention layers to compute keys
        n_units : int
            Number of units in each attention layer
        """
        super(AttentionNetwork, self).__init__()
        self.n_attention_layers = n_layers
        self.n_units = n_units

        fc_attention = []
        for _ in range(self.n_attention_layers):
            att_linear = nn.Linear(n_input_features, self.n_units)
            att_linear.weight.data.normal_(0.0, np.sqrt(1 / np.prod(att_linear.weight.shape)))
            fc_attention.append(att_linear)
            fc_attention.append(nn.SELU())
            n_input_features = self.n_units

        att_linear = nn.Linear(n_input_features, 1)
        att_linear.weight.data.normal_(0.0, np.sqrt(1 / np.prod(att_linear.weight.shape)))
        fc_attention.append(att_linear)
        self.attention_nn = torch.nn.Sequential(*fc_attention)

    def forward(self, inputs):
        """Apply single-head attention network.

        Parameters
        ----------
        inputs: torch.Tensor
            Torch tensor of shape (n_sequences, n_input_features)

        Returns
        ---------
        attention_weights: torch.Tensor
            Attention weights for sequences as tensor of shape (n_sequences, 1)
        """
        attention_weights = self.attention_nn(inputs)
        return attention_weights


class OutputNetwork(nn.Module):
    def __init__(self, n_input_features: int, n_output_features: int = 1, n_layers: int = 1, n_units: int = 32):
        """Output network (`o()` in paper) as fully connected network

        See `deeprc/examples/` for examples.

        Parameters
        ----------
        n_input_features : int
            Number of input features
        n_output_features : int
            Number of output features
        n_layers : int
            Number of layers in output network (in addition to final output layer)
        n_units : int
            Number of units in each attention layer
        """
        super(OutputNetwork, self).__init__()
        self.n_layers = n_layers
        self.n_units = n_units

        output_network = []
        for _ in range(self.n_layers - 1):
            o_linear = nn.Linear(n_input_features, self.n_units)
            o_linear.weight.data.normal_(0.0, np.sqrt(1 / np.prod(o_linear.weight.shape)))
            output_network.append(o_linear)
            output_network.append(nn.SELU())
            n_input_features = self.n_units

        o_linear = nn.Linear(n_input_features, n_output_features)
        o_linear.weight.data.normal_(0.0, np.sqrt(1 / np.prod(o_linear.weight.shape)))
        output_network.append(o_linear)
        self.output_nn = torch.nn.Sequential(*output_network)

    def forward(self, inputs):
        """Apply output network to `inputs`.

        Parameters
        ----------
        inputs: torch.Tensor
            Torch tensor of shape (n_samples, n_input_features).

        Returns
        ---------
        prediction: torch.Tensor
            Prediction as tensor of shape (n_samples, n_output_features).
        """
        predictions = self.output_nn(inputs)
        return predictions