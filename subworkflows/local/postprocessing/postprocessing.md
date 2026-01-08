# POSTPROCESSING Subworkflow

This workflow performs post-alignment processing including BAM sorting, generating allSites, collapsing insertions to sitesFinal, and creating points for clonal expansion analysis.

## Included Modules
- `INDEX_SORT_BAM_local` from `../../../modules/local/index_sort_bam_local/main`
- `BAM_TO_ALLSITES_local` from `../../../modules/local/bam_to_allsites_local/main` 
- `ALLSITES_TO_SITESFINAL_edited_grouping_local` from `../../../modules/local/allsites_to_sitesfinal_local/main`
- `SITESFINAL_TO_POINTS_local` from `../../../modules/local/sitesfinal_to_points_local/main`

The description of the input and output channel structures of each module, along a more detailed description of the workflow, is provided in the specialized Markdown document for each module.

## Workflow: `POSTPROCESSING_wfl`

### Inputs
| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_aligned` | Aligned BAM files from the alignment workflow | Samplename, BAM file |
| `ch_processing_params` | Parameters for post-processing (e.g., filtering thresholds) | minPctIdent, maxAlignStart, maxFragLength
| `ch_refGenome` | Reference genome information for annotation | referenceGenome_name, referenceGenomeFile, referenceKnowngeneFile |

### Main Steps
1. **BAM Sorting**
   - Process: `INDEX_SORT_BAM_local`
   - Input: `ch_aligned`
   - Output channel: `ch_sorted`

   * Sorts BAM files by coordinate and indexes them with samtools.

2. **Prepare BAM input**
   - Merge `ch_sorted` with `ch_processing_params`
   - Output channel: `ch_baminput`
   ```groovy
   ch_baminput = ch_sorted.merge(ch_processing_params)
   ```

3. **Convert BAM to allSites**
   - Process: `BAM_TO_ALLSITES_local`
   - Input: `ch_baminput`
   - Output channel: `ch_allsites`

   * Converts BAM files to allSites format, with one insertion per line.

4. **Collapse allSites to sitesFinal**
   - Process: `ALLSITES_TO_SITESFINAL_edited_grouping_local`
   - Input: `ch_allsites`
   - Output channel: `ch_sitesfinal`

    * Collapses allSites into sitesFinal by grouping duplicated insertions.

5. **Prepare annotation input**
   - Merge `ch_sitesfinal` with `ch_refGenome`
   - Output channel: `ch_annotinput`
   ```groovy
   ch_annotinput = ch_sitesfinal.merge(ch_refGenome)
   ```

6. **Collapse insertions and generate points**
   - Process: `SITESFINAL_TO_POINTS_local`
   - Input: `ch_annotinput`
   - Output channel: `ch_points`

    * Collapses insertions at the same genomic point (same R2 start, same viral insertion) to generate points for clonal expansion analysis.

### Outputs
| Channel | Description | Structure |
|---------|-------------|-----------|
| `ch_points` | Collapsed insertion points, counts, and clonal expansion information | Samplename, plots file, annotation excel

```
