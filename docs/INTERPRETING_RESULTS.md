# Interpreting ZincSight predictions

ZincSight predicts potential transition-metal binding sites from protein structures. Its outputs include predicted metal positions, coordinating residues and site probabilities.

## Metal identity and PyMOL labels

Although trained on zinc-binding sites, ZincSight also detects sites for other transition metals, including iron, manganese, cobalt, copper and nickel, as described in the [published study](https://doi.org/10.1002/pro.70350) (section 2.5, Figure 5 and Discussion).

A predicted site identifies a candidate metal-binding location. ZincSight does not determine which metal species or oxidation state binds that site, or establish that the site binds every metal listed above. Predicted ions labelled Zn in generated PyMOL sessions mark predicted positions; the labels do not establish zinc specificity. Keep generated predictions separate from experimentally resolved metal ions already present in an input structure.

## Site probabilities

The reported probability comes from the zinc-trained scoring model. It is not a set of separately calibrated probabilities for zinc, iron, manganese, cobalt, copper and nickel. Assigning a particular metal requires independent evidence.

When interpreting results with an AI assistant, retain these distinctions and preserve the original output labels and values. A request about iron, for example, should not cause predicted Zn labels to be relabelled as confirmed iron sites.
