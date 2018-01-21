# GeometricFlows

Geometric flows are differential equations that reshape manifolds in order to
simplify or understand their structure. Applications of geometric flows include
finding minimal energy surfaces (mean curvature flow, soap films, etc.), the
proof of the Poincare conjecture (Ricci flow), and machine learning algorithms
The goal of this project is to explore geometric flows using simple examples and
numerical experimentation. The project will include an individual study course
introducing some needed concepts such as partial differential equations,
manifolds, and curvature. Concurrently, the team will work together to develop
numerical algorithms for simulating geometric flows, starting from simple
examples. Once these tools are developed, we will explore various flows by
experimentation and visualization, and in particular we will be interested in
the long time limits of the flow.

---

ProjectTemplate is an R package that helps you organize your statistical
analysis projects. To load this project in R, you'll first need to `setwd()`
into the directory where this README file is located. Then you need to run the
following two lines of R code:

	library('ProjectTemplate')
	load.project()

After you enter the second line of code, you'll see a series of automated
messages as ProjectTemplate goes about doing its work. This work involves:
* Reading in the global configuration file contained in `config`.
* Loading any R packages you listed in he configuration file.
* Reading in any datasets stored in `data` or `cache`.
* Preprocessing your data using the files in the `munge` directory.

Once that's done, you can execute any code you'd like. For every analysis
created, there's a separate file in the `src` directory.
If the files start with the two lines mentioned above:

	library('ProjectTemplate')
	load.project()

You'll have access to all of your data, already fully preprocessed, and
all of the libraries you want to use.

For more details about ProjectTemplate, see http://projecttemplate.net
