# Example of using trelliscope to fit sinusoidal curves to time-series data.

This repo is based on a project that fit cosinor regression models to time-series lipidomics data.  The data used here is fake.

The data is provided in ./data but you can see how it was created in munge_data.R.  There are multiple 'subjects' that are assigned to either the 'day' or 'night' group.  Each sample indicates the group, individual, and time it was taken.  For each subject, they are assigned a sinusoidal response across samples for each feature ('lipids').  The amplitude and phase shift of the sinusoid is forced to be dependent on the group assignment.  Additionally, I induce a random subject effect.

In apply_model.R or Notebook.Rmd, we split the data by each lipid and create a plot across all lipids.  This is done by fitting a cosinor model to the data within each lipid, that essentially tries to align a sinusoid to the data representing each subject for that lipid.  The model fit is then used to draw a plot for each lipid, all of which are arranged in an interactive trelliscope display.
