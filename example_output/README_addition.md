````markdown
#### Input Modes

You can control how sequences are taken into the run CSV with the `--mode` flag:

```bash
--mode {hybrid,pdb_only,seq_only_csv}
````

-----

`pdb_only` (default)

Sequences are extracted **only from the PDB files**.
The run CSV will be updated to reflect what is found in the structures.
Any sequences specified in the CSV are ignored.

-----

`seq_only_csv`

Sequences are taken **only from the CSV**.
Useful if you want full control over sequences (e.g., curated or designed sequences).

-----

`hybrid`

In this mode, by default, sequences are extracted from the PDBs.
However one can also ***overwrite the PDB extracted seqeunces***. To overwrite the PDB extracted sequences with custom sequence, the input CSV must specify both:

  * `target_chains` please specify all target chains (also ones that should not be overwritten)
  * `target_subchain_X_seq` (one column per subchain that one wants to overwrite, e.g., `target_subchain_D_seq`)

Then the CSV-provided sequence will overwrite the extracted version for that chain.

**Additional behaviors in hybrid mode:**

  * **Validation**
    CSV-provided sequences are checked to ensure they contain only valid protein characters: `ACDEFGHIKLMNPQRSTVWY`.
  * **Logging extracted sequences**
    When a CSV-provided sequence overwrites a PDB-extracted one, the extracted sequence is still stored in a dedicated column: `pdb_extracted_trg_subch_{X}_not_used`.
    This ensures transparency: you can always compare the provided sequence against what was found in the structure.

Please see the `example_submission_overwrite` for an example

<!-- end list -->

```
```