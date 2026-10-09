# Upload the tested release

The prepared checkout retains the existing GitHub history and adds one ordinary commit replacing the active workflow. No force push or history deletion is required. Authentication and the final push are performed by the repository owner.

From the prepared checkout on the server:

```bash
cd /absolute/path/github_release_20261009
git remote -v
git log -1 --oneline
git status --short
git push origin HEAD:main
```

If another commit has appeared on GitHub meanwhile, Git rejects the push safely. Fetch and reconcile the changes normally; do not force push. The earlier STAR-based implementation remains accessible through the previous commits.

The code ZIP is convenient for reviewing or uploading the contents through GitHub's file interface after extracting it. Uploading the ZIP alone as a repository file does not replace the pipeline source. The Git bundle can also recreate the prepared repository on another computer; clone the bundle, set origin to `https://github.com/lojour-sketch/inspiired_nf_optimization_v2.git`, then push normally.

Create a GitHub release for the committed version and attach the tested `bushman_runtime.sif` separately, along with its SHA256 file and the code ZIP. Keep the SIF out of the Git source tree. Until that artifact is distributed, the code release requires obtaining the image directly from the owner. No raw FASTQ/BCL/reference sequence or full work directory is included in the code package.

The included reports distinguish the corrected 0% validation from all seven older provisional samples. After a complete corrected RUN809 run, update the biological comparison and stage-runtime evidence; the current release must not claim all-seven corrected equivalence or a complete fresh BCL benchmark.
