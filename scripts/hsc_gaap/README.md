The flow of running `gaap` on HSC S20A data

```mermaid
graph TD;
    Download Images --> Run `gaap.GaapTask` for each patch --> Merge `objectTable` from HSC and Merian, remember to overwrite RA, DEC;
    Merge `objectTable` from HSC and Merian, remember to overwrite RA, DEC --> Concatenate `objectTable` for all patches in a tract;
```