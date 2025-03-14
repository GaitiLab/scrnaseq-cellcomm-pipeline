process testRscript {
    script:
    """
    echo \${PATH}
    /usr/local/bin/_entrypoint.sh Rscript
    """
}
