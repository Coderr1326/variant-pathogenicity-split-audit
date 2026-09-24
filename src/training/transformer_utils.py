"""Shared transformer loading logic for local and Colab execution."""
from torch import nn
from transformers import AutoTokenizer, AutoModel

def load_transformer(model_name: str):
    tokenizer=AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    if model_name == 'zhihan1996/DNABERT-2-117M':
        from transformers import AutoConfig
        cfg=AutoConfig.from_pretrained(model_name, trust_remote_code=True)
        if getattr(cfg, 'pad_token_id', None) is None:
            cfg.pad_token_id=tokenizer.pad_token_id if tokenizer.pad_token_id is not None else 0
        cfg.attention_probs_dropout_prob=0.1
        encoder=AutoModel.from_pretrained(model_name, config=cfg, trust_remote_code=True, low_cpu_mem_usage=False)
    else:
        from transformers import AutoConfig, AutoModelForMaskedLM
        cfg=AutoConfig.from_pretrained(model_name, trust_remote_code=True)
        encoder=AutoModelForMaskedLM.from_pretrained(model_name, config=cfg, trust_remote_code=True, low_cpu_mem_usage=False)
    class SequenceWrapper(nn.Module):
        def __init__(self, model):
            super().__init__(); self.m=model; self.h=nn.Linear(model.config.hidden_size,2)
        def forward(self, **kw):
            kw['output_hidden_states']=True; out=self.m(**kw)
            if hasattr(out,'last_hidden_state'): z=out.last_hidden_state
            elif hasattr(out,'hidden_states') and out.hidden_states is not None: z=out.hidden_states[-1]
            else: z=out[0]
            return self.h(z[:,0])
    return tokenizer, SequenceWrapper(encoder)

def compute_max_length(tokenizer, texts, safety_margin=4, round_to=8):
    """Exact worst-case tokenized length for `texts` under `tokenizer`.

    DNABERT-2's BPE tokenizer produces roughly one token per input character, so a
    character-based max_length (2*window) is a safe, near-tight bound for it. NT's
    tokenizer instead maps every non-overlapping 6-mer to a single token (falling
    back to one token per leftover base), so the same 2*window value pads every
    NT batch to ~5-6x its real token length for no benefit. Pass the few longest
    ref+mut strings by character count (cheap to find, no tokenizer call needed)
    to get the true worst case cheaply instead of reusing DNABERT-2's convention
    or re-tokenizing the whole dataset.
    """
    longest_tokens = max(len(tokenizer(t, add_special_tokens=True)['input_ids']) for t in texts)
    return ((longest_tokens + safety_margin + round_to - 1) // round_to) * round_to
