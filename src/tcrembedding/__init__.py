from .version import __version__

from .ATMTCR.embedding import EmbeddingATMTCR
from .catELMo.embedding import EmbeddingcatELMo
from .clusTCR.embedding import EmbeddingclusTCR
from .DeepRC.embedding import EmbeddingDeepRC
from .DeepTCR.embedding import EmbeddingDeepTCR
from .ERGOII.embedding import EmbeddingERGO
from .ESM.embedding import EmbeddingESM
from .GIANA.embedding import EmbeddingGIANA
from .ImRex.embedding import EmbeddingImRex
from .iSMART.embedding import EmbeddingiSMART
from .LuuEtAl.embedding import EmbeddingLuuEtAl
from .NetTCR2_0.embedding import EmbeddingNetTCR2
from .pMTnet.embedding import EmbeddingpMTnet
from .SETE.embedding import EmbeddingSETE
from .TCRanno.embedding import EmbeddingTCRanno
from .TCRGP.embedding import EmbeddingTCRGP
from .TITAN.embedding import EmbeddingTITAN
from .Word2Vec.embedding import EmbeddingWord2Vec

__all__ = ['EmbeddingATMTCR', 'EmbeddingcatELMo', 'EmbeddingclusTCR', 'EmbeddingDeepRC', 'EmbeddingDeepTCR', 'EmbeddingERGO',
           'EmbeddingESM', 'EmbeddingGIANA', 'EmbeddingImRex', 'EmbeddingiSMART', 'EmbeddingLuuEtAl', 'EmbeddingNetTCR2', 
           'EmbeddingpMTnet','EmbeddingSETE', 'EmbeddingTCRanno', 'EmbeddingTCRGP', 'EmbeddingTITAN', 'EmbeddingWord2Vec']


name = "TCRembedding"