---
layout: post
title: Python, R, Bash in one Jupyter Notebook
date: 2018-12-20 15:26:46 +01:00
categories: Programming
tags: Bash, Jupyter, pipeline, python, R
image: /assets/posts/2018-12-20-python-r-bash-jupyter-notebook/R_Bash_in_Jupyter_notebooks.jpeg
alt: Python, R, Bash in one Jupyter Notebook
excerpt: Combining Python, R, Bash in one Jupyter Notebook makes tracking of the workflow easier, simplifies sharing and makes you more efficient and professional.
---

<h2>Why Jupyter Notebooks?</h2>

<p>If you read <a href="{{ site.baseurl }}/genomic-spark-tutorial/">my Big Data tutorial</a>, you are already familiar with Databricks notebooks. These notebooks allow combining code from many different programming languages (Scala, Python, etc.) in one notebook. I thought it would be great to set up a similar notebook environment locally on my computer to manage my workflows.</p>

<p>I used to log all my work steps in a text editor. It was a simple and reliable approach, but it was not the most efficient one. In particular, it was laborious to copy all commands I execute in the terminal to the text editor. Moreover, such log files were not user-friendly to share with my colleagues.</p>
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-12-20-python-r-bash-jupyter-notebook/log_workflow_in_txt.jpeg" alt="" class="wp-image-1510" />
<figcaption class="aligncenter">Workflow in a txt file.</figcaption>
</figure>

<p>I knew that  a <a aria-label="I knew that  a Jupyter Notebook is the closest I can get to Databriks workflow. But I knew Jupyter Notebooks only as Python Notebooks. However, after a few minutes searching online, I found out that one can easily combine Python, R, Bash in one Jupyter Notebook and execute code in all three languages without leaving a notebook and without manual changing of the Jupyter engine.
 (opens in a new tab)" href="https://jupyter.org/" target="_blank">Jupyter Notebook</a> was the closest to <a href="{{ site.baseurl }}/processing-genomic-data-apache-spark-big-data-tutorial/">Databriks notebook</a> open-source solution that can be implemented on a desktop. But I knew Jupyter Notebooks only as Python Notebooks. However, after a few minutes searching online, I found out that one can easily combine <a aria-label="I knew that  a Jupyter Notebook is the closest to Databriks notebook open-source solution that can be implemented on a desktop. But I knew Jupyter Notebooks only as Python Notebooks. However, after a few minutes searching online, I found out that one can easily combine Python, R, Bash in one Jupyter Notebook and execute code in all three languages within one notebook.
 (opens in a new tab)" href="https://www.python.org/" target="_blank">Python</a>, <a aria-label="I knew that  a Jupyter Notebook is the closest to Databriks notebook open-source solution that can be implemented on a desktop. But I knew Jupyter Notebooks only as Python Notebooks. However, after a few minutes searching online, I found out that one can easily combine Python, R, Bash in one Jupyter Notebook and execute code in all three languages within one notebook.
 (opens in a new tab)" href="https://www.r-project.org/" target="_blank">R</a>, <a aria-label="I knew that  a Jupyter Notebook is the closest to Databriks notebook open-source solution that can be implemented on a desktop. But I knew Jupyter Notebooks only as Python Notebooks. However, after a few minutes searching online, I found out that one can easily combine Python, R, Bash in one Jupyter Notebook and execute code in all three languages within one notebook.
 (opens in a new tab)" href="https://www.gnu.org/software/bash/" target="_blank">Bash</a> in one Jupyter Notebook and execute code in all three languages within one notebook.</p>

<p>I would like to share with you how to set up such a notebook environment.</p>

<h2>R and Bash in Jupyter cells</h2>

<p>Below, I provide instruction for Linux OS, but I believe these commands are similar if not the same as in other operating systems. </p>

<p>First, install <a href="https://jupyter.org/" target="_blank">Jupyter Notebook</a>. You can install it from your package manager by searching for the <code>jupyter</code> or by using <code>pip</code>:</p>

```bash
python -m pip install jupyter
```

<p>Second, install <a href="https://rpy2.bitbucket.io/" target="_blank">an R interface for Jupyter</a> from your package manager by searching for <code>rpy2</code> package or by using <code>pip</code>:</p>

```bash
pip install rpy2
```

<p>Third, load <code>rpy2</code> in your Jupyter Notebook: </p>

```bash
%load_ext rpy2.ipython
```

<p>In my case, it complained that I had to install <code>pandas</code> as a dependency of <code>rpy2.</code> When I installed it, everything worked correctly.</p>

<p><em><strong>Note:</strong> All these packages should be for the same Python version. So, keep track of whether you install these packages for Python 2 or Python 3.</em></p>

<p>After these simple steps, I was able to execute Python, R, Bash in one Jupyter Notebook by indicating R and Bash cells with <code>%%R</code> and <code>%%bash</code>, <a href="https://ipython.readthedocs.io/en/stable/interactive/magics.html" target="_blank">Jupyter magic commands</a>.</p>
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-12-20-python-r-bash-jupyter-notebook/jupyter_notebooks_python-R-bash.jpeg" alt="" class="wp-image-1504" />
<figcaption class="aligncenter">Python, R, Bash in one Jupyter Notebook</figcaption>
</figure>

<p>To test your installation, you can replicate my commands from the image above.</p>

<h2>R and Bash Jupyter kernels</h2>

<p>Basically, if you want to use Jupyter Notebooks primarily for Python code with an option to execute Bash and R, the steps described above are enough. However, you can also go further. If you install R and Bash Jupyter kernels, you will be able to use Jupyter Notebooks as notebooks for pipelines in either of these two languages.</p>

<p>To install <a href="https://github.com/takluyver/bash_kernel" target="_blank">a Jupyter kernel for Bash,</a> execute:</p>

```bash
pip install ipykernel
pip install bash_kernel
python -m bash_kernel.install
```

<p>To install <a aria-label="To install an R kernel for Jupyter, you need to install
 (opens in a new tab)" href="https://github.com/IRkernel/IRkernel" target="_blank">an R kernel for Jupyter</a>, you need to run this code inside R:</p>

```bash
install.packages('IRkernel')
IRkernel::installspec()
```

<p>Now, you should be able to select a particular Jupyter kernel and create a Jupyter Notebook in either Python, Bash or R.</p>
<figure class="caption"><img src="{{ site.baseurl }}/assets/posts/2018-12-20-python-r-bash-jupyter-notebook/jupyter_R_bash_kernels.jpeg" alt="" class="wp-image-1505" />
<figcaption class="aligncenter">Python, Bash, R Jupyter kernels.</figcaption>
</figure>

<p>I usually use R and Bash kernels when I work on exclusive R or Bash pipelines.</p>

<h2>Conclusion</h2>

<p>This way you can use Jupyter Notebooks to log and execute your Python, R, Bash together in one single notebook as well as to create a well-annotated dedicated Python, R, Bash pipelines.</p>

