---
layout: post
title: Snakemake checkpoint tutorial
date: 2019-08-16
categories: Pipelines
tags: Snakemake
image: /assets/posts/2019-08-16-snakemake-checkpoint-tutorial/snakemake-checkpoint-tutorial_thumbnail.jpg
_focus_key_word: Snakemake checkpoint tutorial
excerpt: "This Snakemake checkpoint tutorial shows how to run Snakemake when the number of outputs is dynamic, e.g. file names are unknown until the rule is executed."
---

If you want to use Snakemake to run some programs that output an unknown number of files, you need to tell Snakemake about that. If you use Snakemake 4, you can do that by marking the output with `dynamic()`. If you upgraded to Snakemake 5, you better use `checkpoint`. Using `dynamic()` will work in Snakemake 5, but you will see a message saying that dynamic output is deprecated and will be fully replaced by checkpoints in Snakemake 6.

This post shows how to use both `dynamic()` and `checkpoint`.

You probably better focus on `checkpoint` because this is a more up-to-date solution. But `checkpoint` may not work correctly sometimes. For example, I tested it with the GATK `IntervalListTools` and [it did not work correctly](https://stackoverflow.com/questions/57432036/snakemake-checkpoint-exited-with-non-zero-exit-code){:target="_blank"}, while `dynamic()` worked fine with the exactly same command. Thus, knowing both approaches can be helpful.

## Checkpoint

[Checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution){:target="_blank"} function was introduced in Snakemake 5 and it will completely replace `dynamic()` in Snakemake 6. So, if you have not tried it, it is time to learn it.

Here is a dummy code that shows how `checkpoint` works:
```python
rule final_output:
    input:
        'scatter_copy_head_collect/all.txt'

# generate random number of files
checkpoint scatter:
    output:
        directory('scatter')
    shell:
        '''
        N=$(( $RANDOM % 10))
        for j in $(seq 1 $N); do echo -n $j > scatter/$j.txt; done
        '''

# process these unknown number of files
rule scatter_copy:
    output:
        txt = 'scatter_copy/{i}_copy.txt',
    input:
        txt = 'scatter/{i}.txt',
    shell:
        '''
        cp -f {input.txt} {output.txt}
        echo -n "_copy" >> {output.txt}
        '''
# process scatter_copy output
rule scatter_copy_head:
    output:
        txt = 'scatter_copy_head/{i}_head.txt',
    input:
        txt = 'scatter_copy/{i}_copy.txt',
    shell:
        '''
        cp -f {input.txt} {output.txt}
        echo "_head" >> {output.txt}
        '''

# collect the results of processing unknown number of files
# and merge them together into one file:

def aggregate_input(wildcards):
    '''
    aggregate the file names of the random number of files
    generated at the scatter step
    '''
    checkpoint_output = checkpoints.scatter.get(**wildcards).output[0]
    return expand('scatter_copy_head/{i}_head.txt',
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.txt')).i)

rule scatter_copy_head_collect:
    output:
        combined = 'scatter_copy_head_collect/all.txt',
    input:
        aggregate_input
    shell:
        '''
        cat {input} > {output.combined}
        '''
```

Explore the outputs, to understand how this pipeline works:

![unknown output files of snakemake with checkpoint rule](/assets/posts/2019-08-16-snakemake-checkpoint-tutorial/snakemake-checkpoint-output.jpg)

## Dynamic

[Dynamic output](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#dynamic-files){:target="_blank"} is outdated approach but it seems to be more stable and reliable at the moment. So, if you experience some problems with `checkpoint`, in most cases, you can write the same pipeline with `dynamic()`.

This is the same pipeline as above but it utilizes `dynamic()` instead of `checkpoint`:

```python
rule final_output:
    input:
        'scatter_copy_head_collect/all.txt'

# this was a checkpoint step above:
rule scatter:
    output:
        dynamic('scatter/{i}.txt')
    shell:
        '''
        N=$(( $RANDOM % 10))
        for j in $(seq 1 $N); do echo -n $j > scatter/$j.txt; done
        '''

# this rule is not different from checkpoint
rule scatter_copy:
    output:
        txt = 'scatter_copy/{i}_copy.txt',
    input:
        txt = 'scatter/{i}.txt',
    shell:
        '''
        cp -f {input.txt} {output.txt}
        echo -n "_copy" >> {output.txt}
        '''

# this rule is not different from checkpoint either:
rule scatter_copy_head:
    output:
        txt = 'scatter_copy_head/{i}_head.txt',
    input:
        txt = 'scatter_copy/{i}_copy.txt',
    shell:
        '''
        cp -f {input.txt} {output.txt}
        echo "_head" >> {output.txt}
        '''

# to collect all files, you need to tell Snakemake that input is dynamic:
rule scatter_copy_head_collect:
    output:
        combined = 'scatter_copy_head_collect/all.txt',
    input:
        indivfiles = dynamic('scatter_copy_head/{i}_head.txt')
    params:
        gathered = lambda wildcards, input: ' '.join(input.indivfiles)
    shell:
        '''
        cat {params.gathered} > {output.combined}
        '''
```

## Final thoughts

Checkpoints are claimed to be more powerful that `dynamic()` by the Snakemake developers. I believe they are right but my impression is that `dynamic()` is easier to use. Maybe I have not fully comprehended `checkpoint` yet.

Besides, as I mentioned above I was not able to make it work with GATK. So, I will try to use `checkpoint` but I may also step back to `dynamic()` too. 

Finally, I would like to acknowledge this [Stackoverflow answer](https://stackoverflow.com/a/56451259/2317701){:target="_blank"} that inspired me to write this tutorial. 
