# ali_imaging

This repository contains code supporting the manuscript "Live imaging of airway epithelium reveals that mucociliary clearance modulates SARS-CoV-2 spread." Code used in data processing, data analysis, and figure generation is provided. 
* R project & data files used for figure generation and data analysis are provided in 'figures'.
 * * 'data' contains data that were small enough to be uploaded here, not in the dryad repository.
   * 'exptlog_iii.csv' is annotations & some measurements (piv, thresholding for gfp and spy650-tubulin, qPCR at 120 hpi, etc) for each live imaged culture.
   * 'hunting_ii.csv' is single cell tracking notes. 'venus_omicron_survival.csv' is as well, for cultures infected with those reporters.
   * Anything with 'qpcr' in the name is N concentration in apical rinsate.
   * '240702_goodspots.csv' is collated spot detection data from live imaged cultures.
   * 'raw_qpcrdata' contains unprocessed data files for each run.
   * Other data are called as used in the .R files.
* Conda environments used for image analysis and processing are provided in 'environments'.
* 'ciliary_beat_frequency' contains code used for analyzing ciliary motion movies.
  * 'denoise.ipynb' implements cellpose 3 denoising of ciliary motion movies.
  * 'beatfreq.ipynb' processes ciliary motion movies (denoised or raw), calculating beat frequency of each pixel, collating with reference image data if available, and concatenating output .csv files into more manageable per-experiment .csv files.
  * 'pts_exploration.R' is useful for determining exactly which points correspond to which culture in a given experiment, based on .pts files saved during image acquisition.
* 'whole_video_processing' contains code used for processing whole culture movies.
  * 'stitchem5.ijm' stitches raw .dv movie tilescans. It is an imageJ macro. Parameters were optimized for this specific project.
  * 'cv2_circles.ipynb' automatically zeroes out non-culture area pixels.
  * 'piv.ipynb' performs particle image velocitometry on the first channel of fully processed whole culture movies. Used for quantifying cell migration.
  * 'thresholding.ipynb' was used for signal intensity quantification of whole culture movies.
  * 'spots_file_240725.r' was used to collate Trackmate spot detection data derived from fully processed whole culture movies.

 # Contact
 If you have any questions, reach out to the corresponding author Thomas Hope at thope@northwestern.edu or the author of this repository Mark Becker at mbecker@u.northwestern.edu.
