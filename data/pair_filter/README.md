# miRNA–TF pair filter: inputs and pair-level output

These files are the inputs and output of the miRNA–TF pair filter of the first analysis round (Methods 2.1 of the
paper; Supplementary Note S1). They are deposited as they were produced, except that the workbook's document
properties now name the authors of the paper and Excel's record of the folder it was last saved in was removed;
no cell of the workbook was changed. `universe_256_genes.txt` was added later: it lists the genes identified as the
filter's universe (below).

The filter did not define the network, and no other result depends on these files:
- the network itself is in `data/network/`;
- every FFL analysed in the paper is formed by edges of that network, not by pairs of this filter;
- mapped to the network's nodes, 11 of the network's 223 miRNAs and 42 of its 157 TF-typed nodes occur only in pairs the
  filter rejected;
- the step that assembled the network's source FFL networks was not recorded.

| File | Content |
|---|---|
| `Hypergeometric_Test_2.xlsx` | Pair-level output, one sheet. The left block has, for each tested row, the miRNA, the TF, `Nmir`, `Ntf`, `i` (the number of genes targeted by both) and the *P* value. The right block has the same *P* values ranked, with Benjamini–Hochberg critical values, adjusted *P* values and the retention call ("Significant using an FDR of 0.05?"). |
| `value.txt` | One line per tested row: miRNA, TF, `Nmir`, `Ntf`, `i`. This is the input block of the workbook. |
| `mirna-tf.txt` | The tested miRNA–TF pairs: TransmiR human TF→miRNA records. |
| `common-i.txt` | One line per tested row: miRNA, TF, `i`. |
| `Nmir.txt` | Target count per miRNA. |
| `Ntf.txt` | Target count per TF; its 795 regulators are those of TRRUST v2 (human). |
| `mirna-gene.txt` | miRNA–target pairs (source not recorded). |
| `tf-gene.txt` | TF–target pairs: TRRUST v2 (human). |
| `universe_256_genes.txt` | The 256 genes targeted in both `mirna-gene.txt` and `tf-gene.txt`, taken to be the filter's universe (added; see below). |

## How the *P* value was computed

The *P* value in the workbook is the lower tail of the hypergeometric distribution, *P*(*X* ≤ *i*), with a
universe of *N* = 256 genes. A pair was retained when its Benjamini–Hochberg-adjusted *P*, computed over all
2,676 rows, was below 0.05.

The test therefore asks whether a pair shares fewer targets than expected, not more:
- The 2,676 rows hold 2,253 distinct tested pairs, and the 1,865 retained rows hold 1,566 distinct pairs.
- Of the 1,566 retained pairs, 1,547 share no target, and none shares more targets than expected.
- None of the 2,111 distinct pairs with a defined *P* value reaches *P* < 0.05 in the upper tail.

## Known features of the files

**`value.txt` cannot be rebuilt from the other files, because its inputs are inconsistent:**
- `Nmir` is 85 in every row.
- 142 distinct pairs (181 rows) have `Ntf` values above 256, so the workbook shows "#NUM!" for them. Eighteen of the
  20 regulator entries involved are absent from `Ntf.txt` and carry the count of *NFKB1* (303) or *RELA* (301).
- `common-i.txt` differs from `value.txt` in 3 rows; the workbook used `value.txt`.

**The 2,676 rows include repeats and invalid rows.** They contain 423 rows that repeat a pair and the 181 "#NUM!"
rows. All 2,676 rows entered the Benjamini–Hochberg correction.

**The universe.** No file of the first round names the 256 genes. Exactly 256 genes are targeted in both
`mirna-gene.txt` and `tf-gene.txt`, so these are taken to be the universe; the identification rests on this count.
They are listed in `universe_256_genes.txt`.

**Entries in the retained pairs:**
- 344 miRNA entries, which are precursor or family names;
- 236 regulator entries, of which 131 are among the 795 TRRUST regulators of `Ntf.txt`. The others entered through
  the miRNA–TF pair list and include signalling proteins and non-gene labels, for example TGFB1, IL6 and LPS.

**Spreadsheet artefacts:**
- DEC1 appears as a date or as the serial number 37226.
- Four TF–target entries in `tf-gene.txt` have their target symbol converted to a date ('07-Sep' once, '01-Dec' three
  times).
- One miRNA name in `mirna-gene.txt` carries a duplicated leading "h" (hhsa-miR-211-5p).

## MD5 checksums

| File | MD5 |
|---|---|
| `Hypergeometric_Test_2.xlsx` | `261b35d2cb6638a84a29e6f3c834eb3a` |
| `Nmir.txt` | `5d2ba5c05151761a643001521acc2503` |
| `Ntf.txt` | `d6212b51e6b54e966bc905a469b78980` |
| `common-i.txt` | `361c2c8f7e230de356ed8f2d28c831aa` |
| `mirna-gene.txt` | `ec4c13aa693f8078ff5dfe06171d3927` |
| `mirna-tf.txt` | `e024c284a4b7d67d6a9b9c21ac07d7ba` |
| `tf-gene.txt` | `28d47380778d404e083f90eab6efc59e` |
| `value.txt` | `54d4311c318713eba806a35158057908` |
| `universe_256_genes.txt` | `dc7b59f4ac71c314594b5d95e67704bf` |
