# GitHub -> Zenodo release checklist

Before making the repository public and citing it in the manuscript:

1. **Check privacy**
   - Run `git status` and inspect every staged file.
   - Confirm no DICOM/NIfTI/MGZ/MAT/CSV/log/FreeSurfer subject folders are present.
   - Search for names, hospital identifiers, MRNs, dates of birth and local patient IDs.

2. **Check local paths**
   - Search for `/Users/`, personal home-directory names and private server paths.

3. **Freeze the environment**
   - Save the Python environment used for the paper, for example:
     `conda env export --from-history > environment_from_history.yml`
   - Save pip-resolved versions if relevant:
     `pip freeze > pip-freeze.txt`
   - Record MATLAB, SPM, FreeSurfer and `stroke_segmentor` versions in the GitHub release notes.

4. **Run code checks**
   - `python -m py_compile` on all Python scripts.
   - `zsh -n` on shell scripts.
   - Run the final analysis on the frozen private dataset and verify that key tables/figures reproduce the submitted results.

5. **Choose a software license**
   - Add a license before public release. MIT is a common permissive option for academic analysis code, but the choice is the author's.

6. **Update citation metadata**
   - Confirm author name(s) in `CITATION.cff`.
   - Add ORCID if desired.
   - Set the release version and date.

7. **Create GitHub release**
   - Suggested first manuscript release: `v1.0.0`.
   - Use a release title that matches the manuscript analysis version.

8. **Archive on Zenodo**
   - Connect the GitHub repository to Zenodo.
   - Archive the exact `v1.0.0` release.
   - Copy the resulting DOI.

9. **Insert DOI into the paper**
   - Update `CITATION.cff` with the DOI.
   - Add the DOI to the manuscript Methods/Data Availability Statement.

10. **Do not rewrite the archived release**
    - Any later correction should become `v1.0.1`, `v1.1.0`, etc.
