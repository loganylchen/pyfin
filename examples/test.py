from fin.io.io_pod5 import Pod5Reader
from fin._f5c import detect_events, eventalign, profile_hmm_eventalign, is_available
import pandas as pd

pod5_reader = Pod5Reader("../../testdata/test.pod5")
pod5_reader.open()
read = pod5_reader.get_read("4816bde4-c87f-4eaf-96b9-2c8d3ce8d5d1")
df = pd.DataFrame(detect_events(read.signal_pa))
seq = "TCAACCGGGTTTTAC"
result = profile_hmm_eventalign(read.signal_pa, seq[::-1], 5)
eventalign_df = pd.DataFrame(result["alignment"]).loc[
    :,
    [
        "ref_position",
        "ref_kmer",
        "event_idx",
        "model_kmer",
        "event_mean",
        "event_duration",
        "model_mean",
        "scaled_model_mean",
    ],
]
df.loc[(df["mean"] > 79.85) & (df["mean"] < 79.86)]
eventalign_df.loc[(eventalign_df["event_mean"] > 79.85) & (eventalign_df["event_mean"] < 79.86)]
