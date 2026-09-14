# dang_signal_mod

Acts as the component factory/dispatch layer.
Its initialization routines allocate component lists and instantiate concrete component types from parsed `fg_type` settings.

In practice, this is where the configured foreground model becomes actual objects.
It wires user-defined component choices to the corresponding module implementations before fitting starts.

At runtime it builds the model list so total signal is assembled as

$$
\hat d = \sum_c \hat d_c,
$$

with each \(\hat d_c\) coming from the component subtype selected by `fg_type`.
