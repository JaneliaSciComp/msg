#!/usr/bin/env Rscript

## [[file:~/src/dedR/ded.R.org::*Matrix%20operations][block-1]]
read.matrix <- function(file, colClasses="numeric", row.names, col.names, ...) {
    x <- as.matrix(read.table(file=file, colClasses=colClasses, row.names=row.names, col.names=col.names, ...))
    if(missing(col.names)) colnames(x) <- NULL
    if(missing(row.names)) rownames(x) <- NULL
    x
}

read.matrix.fast <- function(file, nrows, ncols, what=numeric()) {
    if(is.list(what)) ncols <- sum(!sapply(what, is.null))
    else what <- rep(list(what), ncols)
    ll <- scan(file, what=what, nlines=nrows, quiet=TRUE)
    array(unlist(ll), dim=c(nrows, ncols))
}

alternate <- function(x, y, margin) {
    stopifnot(dim(x) == dim(y), length(x) == length(y))
    if(is.null(dim(x))) c(rbind(x,y))
    else if(length(dim(x) == 2)) {
        stopifnot(margin %in% 1:2)
        if(margin == 2) { x <- t(x) ; y <- t(y) }
        nr <- dim(x)[1]
        nc <- dim(x)[2]
        ans <- t(structure(c(aperm(abind(x, y, along=3), perm=c(2,3,1))), dim=c(nc, 2*nr)))
        if(margin == 2) t(ans) else ans
    }
    else stop("What does alternate mean for arrays of dimension > 2 ?")
}

sum.matrices <- function(...) {
    ## This should be re-written using Reduce
    xx <- list(...)
    if(is.list(xx[[1]])) xx <- xx[[1]]
    dims <- sapply(xx, dim)
    stopifnot(all(dims == dims[1:2]))
    n <- length(xx)
    if(n > 1) for(i in 2:n) xx[[1]] <- xx[[1]] + xx[[i]]
    xx[[1]]
}

sum.adjacent.columns <- function(x) {
    stopifnot((twon <- ncol(x)) %% 2 == 0)
    n <- twon / 2
    xo <- x[,seq.int(1, twon-1, 2)]
    xe <- x[,seq.int(2, twon,   2)]
    matrix(xo + xe, nrow(x), n)
    ## S <- matrix(c(rep(rep(c(1,0), c(2,twon)), n-1), c(1,1)), nrow=twon, ncol=n)
    ## x %*% S
}

sum.adjacent.rows <- function(x) t(sum.adjacent.columns(t(x)))

group.sums <- function(x, labels, na.rm=FALSE, average=FALSE) {
    labels <- factor(labels)
    stopifnot((n <- length(labels)) == ncol(x))

    group.sums.with.nas <- function(x, labels, average=FALSE) {
        K <- nlevels(labels)
        p <- array(dim=c(nrow(x), K))
        rownames(p) <- rownames(x)
        cn <- colnames(p) <- levels(labels)
        f <- if(average) rowMeans else rowSums
        for(k in 1:K)
            p[,k] <- f(x[,labels == cn[k], drop=FALSE], na.rm=TRUE)
        p
    }

    if(na.rm)
        return(group.sums.with.nas(x, labels, average))

    y <- sapply(levels(labels), function(k) labels == k)
    ans <- x %*% y
    if(average) {
        n <- apply(y, 2, sum)
        ans <- ans %*% diag(1/n)
    }
    dimnames(ans) <- list(rownames(x), levels(labels))
    ans
}

group.means <- function(x, labels, na.rm=FALSE)
    group.sums(x, labels, na.rm=na.rm, average=TRUE)
## block-1 ends here

## [[file:~/src/dedR/ded.R.org::*Matching%20combinations%20labelling][block-2]]
compare.groups <- function(f1, f2, quiet=TRUE) {
    stopifnot((n <- length(f1)) == length(f2))
    cmp <- function(i) {
        if(!quiet && i %% 1000 == 0) cat(i, "\r")
        dist(rbind(f1 == f1[i], f2 == f2[i]), method="binary")
    }
    ans <- sapply(1:n, cmp)
    if(!quiet) cat("\n")
    ans
}

optimize.permn <- function(fun, vec, nreps=n^3, quiet=TRUE) {
    ## Attempt to find minimum of fun(perm) over all permutations of
    ## vec, using nreps independent repeats of a hill-climbing
    ## algorithm. At each step of the hill-climbing algorithm, the
    ## current permutation v is compared with all permutations that
    ## can be obtained by swapping two elements of v. If the objective
    ## function takes a lower value at one of those neighbouring
    ## points, then the algorithm moves to the point with the lowest
    ## such value and continues; if not then the algorithm terminates.

    n <- length(vec)
    stopifnot(n > 1)
    seq1n <- 1:n
    climb <- function() {
        swap <- function(els) { v[els] <- v[rev(els)] ; v }
        v <- sample(vec)
        e <- fun(v)
        if(!quiet) cat(".")
        repeat {
            vv <- utils::combn(seq1n, 2, fun=swap)
            ee <- apply(vv, 2, fun)
            wmin <- which.min(ee)
            enew <- ee[wmin]
            if(enew >= e) break
            e <- enew
            v <- vv[,wmin]
        }
        attr(e, "perm") <- v
        e
    }
    ee <- replicate(nreps, climb(), simplify=FALSE)
    if(!quiet) cat("\n")
    wmin <- which.min(unlist(ee))
    list(minimum=attr(ee[[wmin]], "perm"), objective=c(ee[[wmin]]))
}

## my solution to an r-help question
occur.matrices <- function(df) {
    levels <- levels(unlist(df))
    ans <- lapply(levels, function(level) crossprod(df == level))
    structure(ans, names=levels)
}

matched.sample <- function(template, pool) {
    ## 'template' and 'pool' are (at least conceptually) factors.
    ## For each level of 'template', there are at least as many elements of 'pool' with the same value.
    ## Generate indices specifying a random sample 'x' from pool, that matches 'template',
    ## i.e. such that all(table(x) == table(template))

    tab <- table(template)
    tab <- tab[tab > 0]
    levs <- names(tab)
    msam <- unlist(lapply(levs, function(lev) sample(which(pool == lev), size=tab[lev], replace=FALSE)))
    structure(msam, names = rep(levs, tab))
}

intersect.arrays <- function(xx, along=2) {
    dims <- lapply(xx, dim)
    stopifnot(sapply(dims, length) == 2)
    stopifnot(along %in% 1:2)

    if(1 %in% along) rnn <- rep(list(Reduce(intersect, lapply(xx, rownames))), length(xx))
    else rnn <- lapply(xx, function(x) seq_len(nrow(x)))

    if(2 %in% along) cnn <- rep(list(Reduce(intersect, lapply(xx, colnames))), length(xx))
    else cnn <- lapply(xx, function(x) seq_len(ncol(x)))

    mapply(function(x,i,j) x[i,j,drop=FALSE], xx, rnn, cnn, SIMPLIFY=FALSE)
}
## block-2 ends here

## [[file:~/src/dedR/ded.R.org::*Distance%20matrix%20S3%20class][block-3]]
"[.dist" <- function(d, i, j) {
    n <- attr(d, "Size")
    if(missing(i)) i <- 1:n
    if(missing(j)) j <- 1:n
    if(is.character(i) || is.character(j)) {
        stopifnot(!is.null(labs <- attr(d, "Labels")))
        i <- if(is.character(i)) match(i, labs) else i
        j <- if(is.character(j)) match(j, labs) else j
    }
    stopifnot(i <= n, j <= n)
    ans <- as.vector(d)[dist.index.vec(i, j, n)]
    ans[is.na(ans)] <- 0
    matrix(ans, nrow=length(i), byrow=TRUE)
}

"[<-.dist" <- function(d, i, j, value) {
    n <- attr(d, "Size")
    if(missing(i)) i <- 1:n
    if(missing(j)) j <- 1:n
    klass <- class(d)
    d <- unclass(d)
    d[dist.index.vec(i, j, n)] <- value
    class(d) <- klass
    d
}

dist.index.vec <- function(i, j, n) {
    ii <- rep(i, each =length(j))
    jj <- rep(j, length(i))
    sapply(seq_along(ii), function(k) dist.index(ii[k], jj[k], n))
}

dist.index <- function(i, j, n) {
    ## Indexing for a "dist" object; see ?dist

    ## I think it follows from considering the values occupying the
    ## lower-triangle of a square matrix with columns labelled i and
    ## rows labelled j. Entry (j,i) is (j-i) values into the ith
    ## column. The number of entries in the preceding (i-1) columns is
    ## choose(n,2) - choose(n-i+1,2).
    stopifnot(length(i) == 1, length(j) == 1)
    if(i == j) return(NA)
    else if(i > j) {
        tmp <- i ; i <- j ; j <- tmp
    }
    stopifnot(j <= n)
    n*(i-1) - i*(i-1)/2 + j-i
}

dist.index.inv <- function(idx, n) {
    ## See explanation in dist.index()

    ## Get indices of topleft values in triangles counting from bottom right
    triangles.rev <- (1:n)*(0:(n-1))/2
    ## Indexing starting from bottom right
    triangles <- rev(triangles.rev[n] - triangles.rev + 1)
    stopifnot(idx <= triangles[n-1])
    ## Which column is the index in?
    i <- max(which(triangles <= idx))
    ## Get j index (how far down column)
    j <- i + idx - triangles[i] + 1
    c(i, j)
}

pairwise.match <- function(d, f, t, nmax, quiet=TRUE) {
    n <- attr(d, "Size")
    ans <- rep(NA, n)
    stopifnot(xor(missing(t),missing(nmax)))
    
    ## Set intragroup distances to NA
    f <- factor(f)
    stopifnot(nlevels(f) == 2)
    klass <- class(d)
    d <- unclass(d)
    for(lev in levels(f)) {
        ii <- which(f == lev)
        d[dist.index.vec(ii, ii, n)] <- NA
    }
    class(d) <- klass
    
    ## Successively pair closest intergroup pairs.
    pair <- 1
    numleft <- table(f)
    while(all(numleft > 0)) {
        minidx <- which.min(d)
        if(!missing(t) && as.vector(d)[minidx] > t) break
        else if(pair == maxpairs) break
        minij <- dist.index.inv(minidx, n)
        ans[minij] <- pair
        d[minij,] <- NA
        pair <- pair + 1
        numleft <- numleft - 1
        #if(!quiet && pair %% 100 == 0) cat(pair, "\r")
        if(!quiet) cat(pair, "\r")
    }
    if(!quiet) cat("\n")
    factor(ans)
}

distance.matrix.summary <- function(dmat, labs, fun=mean, ...) {
    stopifnot(length(labs) == dim(dmat))
    ltri <- lower.tri(dmat, diag=FALSE)

    rlabs <- as.character(labs[row(dmat)][ltri])
    clabs <- as.character(labs[col(dmat)][ltri])

    ## ilabs <- interaction(rlabs, clabs, drop=TRUE)
    ilabs <- factor(paste(pmin(rlabs, clabs), pmax(rlabs, clabs), sep="."))

    dmat <- dmat[ltri]

    vals <- split(dmat, ilabs)
    sapply(vals, fun, ...)
}
## block-3 ends here

## [[file:~/src/dedR/ded.R.org::*Misc][block-4]]
approx.multi <- function(x, y, xout) {
    ## approx when there may be multiple solutions for xout
    ## Josh and Dan
    stopifnot(length(xout) == 1)
    ii <- which(diff(x > xout) != 0)
    y <- sapply(ii, function(i){approx(x[i:(i+1)], y[i:(i+1)], xout=xout)$y})
    list(x=xout, y=y)
}

factor.contiguous <- function(x) {
    ## Assign integer codes to stretches of contiguous positions
    if(FALSE) {
        ## should be able to do it with
        ## rle(diff(x) != 1)
    }
    else {
        n <- length(x)
        fac <- rep(NA, n)
        if(n == 0) return(fac)
        k <- 1
        fac[1] <- k
        if(n == 1) return(fac)
        for(i in 2:n)
            fac[i] <- if(abs(x[i] - x[i-1]) <= 5) k else (k <- k+1) ### arbitrary distance, needed because of indels
            #fac[i] <- if(x[i] - x[i-1] == 1) k else (k <- k+1)
    }
    factor(fac)
}

more <- function(x, chunk=58) {
    if(is.null(dim(x)))
        dim(x) <- c(length(x), 1)
    n <- nrow(x)
    i <- 1
    while(i < n) {
        print(x[i:min(i+chunk,n),,drop=FALSE])
        readline("")
        i <- i + chunk
    }
}

ellipse.contains <- function(xy, rx=1, ry=a, mid=c(0,0), angle=0) {
    ## Return logical vector ans: ans[i] is TRUE iff xy[i,] is inside
    ## an ellipse with angle theta, semimajor axis a, semiminor axis
    ## b, centred at centre.
    xy <- t(t(xy[,1:2]) - mid)
    if(angle != 0) {
        require(shape)
        xy <- shape::rotatexy(xy, angle=-angle, mid=c(0,0))
    }
    colSums( (t(xy)/c(rx,ry))^2 ) < 1
}

mycurve <- function(f, interval, add=FALSE, col=1, ...) {
    x <- seq(min(interval), max(interval), length.out=100)
    y <- sapply(x, f, ...)
    f <- if(add) lines else plot
    f(x, y, type="l", col=col)
}

tooclose <- function(map, d) {
    ## a bit like duplicated()
    sep <- diff(map)
    if(!all(sep > 0)) stop("map must be increasing")
    c(FALSE, sep <= d)
}

all.same <- function(x) all(duplicated(x)[-1])

coordinates <- function(cells, mat) {
    rows <- row(mat)
    cols <- col(mat)
    rnames <- rownames(mat)
    cnames <- colnames(mat)
    if(any(sapply(dimnames(mat), is.null)))
        coord <- function(i) c(rows[i], cols[i])
    else coord <- function(i) c(rnames[rows[i]], cnames[cols[i]])
    t(sapply(cells, coord))
}


write.tabular <- function(x, file=stdout(), ...)
    write.table(x, file=file, row.names=FALSE, col.names=FALSE, quote=FALSE, ...)

lower.diag.proportions <- function(chunks, nrows=16179) {
    if(sum(chunks) != nrows) stop("sum(chunks) falls ", nrows - sum(chunks), " short of ", nrows)
    
    f <- function(n) n*(n+1)/2
    s <- diff(c(0,f(cumsum(chunks))))
    s / sum(s)
}

logmsg <- function(...) {
    cat(date(), "\t", ..., "\n")
}

pad <- pad.chrom <- function(i, width=2) {
    stopifnot(width == 2, i >= 0)
    if(i < 10) paste("0", i, sep="") else as.character(i)
    
}

scan.column <- function(j, file, what=numeric(), ...) {
    command <- paste("awk '{print $", j, "}' ", file, sep="")
    scan(pipe(command), what=what, ...)
}

transpose.file <- function(infile, outfile, ncols, nrows, chunk.size=ncols, what=integer()) {
    stopifnot(!file.exists(outfile))
    if(missing(nrows) || missing(ncols)) {
        tmp <- wc(infile)
        nrows <- tmp["lines"]
        ncols <- tmp["words"] / tmp["lines"]
        stopifnot(floor(ncols) == ncols)
    }
    
    f <- file(outfile, open="wt")
    last <- 0
    while(last < ncols) {
        chunk <- min(chunk.size, ncols - last)
        cat("reading", chunk, "columns:\t", last+1, "to", last+chunk.size, "...")
        x <- scan(infile, what=rep(c(list(NULL), list(what), list(NULL)), c(last, chunk, ncols-chunk-last)), nlines=nrows)
        cat("\nappending columns as rows...")
        lapply(x[!sapply(x, is.null)], function(col) { cat(col, file=f, append=TRUE) ; cat("\n", file=f, append=TRUE) } )
        cat("\n")
        rm(x) ; gc()
        last <- last + chunk.size
    }
    close(f)
}

decode <- function(obj) UseMethod("decode")

decode.factor <- function(f) levels(f)[f]

time.cmp <- function(expr1, expr2, n=10) {
    t <- array(dim=c(n,6))
    for(i in 1:n)
        t[i,] <- c(system.time(parse(text=expr1))[1:3],system.time(parse(text=expr2))[1:3])
    
    ## t(sapply(1:n, function(dummy) c(system.time(expr1)[1 + dummy - dummy],system.time(expr2)[2 + dummy - dummy])))
    t
}

is.even <- function(i) (i %% 2) == 0 ##2 * floor(i/2) == i
## block-4 ends here

## [[file:~/src/dedR/ded.R.org::*Programming][block-5]]
er <- function(on=TRUE) options(error=if(on) recover else NULL)

getopts <- function(argv=NULL) {
    ## for argument-processing with Rscript / littler
    ##
    ## return a list with elements whose names correspond to the flags provided,
    ## and whose values contain the value provided; NULL if no value given.

    if(is.null(argv)) argv <- commandArgs(trailing=TRUE)
    if(!(argc <- length(argv))) return(NULL)
    flagsi <- grep("^-[a-zA-Z]", argv)
    nopts <- length(flagsi)
    valuesi <- setdiff(1:argc, flagsi)
    nvalues <- length(valuesi)
    d <- diff(flagsi)
    stopifnot(d %in% 1:2) ## values, if any, must be a single 'word' (i.e. quote vector values)
    gotval <- d == 2
    ngotval <- sum(gotval)
    if(nvalues == ngotval+1) d <- c(d, 2) ## last flag is followed by a value
    else {
        stopifnot(nvalues == ngotval)
        d <- c(d, 1) ## last flag has no value
    }
    opts <- vector(mode="list", length=nopts)
    names(opts) <- sub("-", "", argv[flagsi])
    opts[d == 2] <- argv[valuesi]
    opts
}

vectorise <- function(f, ...) function(x) sapply(x, f, ...)

get.named.dims <- function(array) structure(dim(array), names=names(dimnames(array)))
## block-5 ends here

## [[file:~/src/dedR/ded.R.org::*R%20environment][block-6]]
lsnof <- function(env=globalenv()) {
    obj.names <- ls(envir=env)
    is.fn <- sapply(obj.names, function(x) is.function(eval(as.name(x), envir=env)))
    obj.names[!is.fn]
}

rmnof <- function() rm(list=lsnof(env=globalenv()))

ls.object.sizes <- function(env=globalenv()) {
    nof <- lsnof(env=env)
    if(length(non.functions) > 0) {
        bytes <- sort(sapply(nof, function(s) object.size(get(s, env=env))), decreasing=TRUE)
        bytes * 2^(-20)
    }
    else nof
}

read.object <- function(f) {
    load(f)
    rm("f")
    eval(parse(text=ls()))
}
## block-6 ends here

## [[file:~/src/dedR/ded.R.org::*System%20tools][block-7]]
wc <- function(f, lines.only=FALSE, gz=FALSE) {
    if(lines.only) {
        call <- if(gz) paste("zcat", f, "| wc -l") else paste("wc -l <", f)
        structure(as.integer(system(call, intern=TRUE)), names="lines")
    }
    else {
        call <- if(gz) paste("zcat", f, "| wc") else paste("wc <", f)
        x <- as.integer(strsplit(gsub(" +", " ", system(call, intern=TRUE)), split=" ")[[1]][-1])
        structure(x, names=c("lines","words","bytes"))
    }
}

battery.info <- function() {
    file <- "/proc/acpi/battery/CMB1/state"
    if(!file.exists(file)) return(Inf)
    txt <- system(paste("cat", file), intern=TRUE)
    state <- grep("charging state", txt, value=TRUE)
    state <- strsplit(state, split = "\\s+")[[1]][3]
    if(state == "charging") hours.remaining <- Inf
    else if(state == "discharging") {
        cap <- grep("remaining capacity", txt, value=TRUE)
        mAh <- as.numeric(strsplit(cap, split="\\s+")[[1]][3])
        rate <- grep("present rate", txt, value=TRUE)
        mA <- as.numeric(strsplit(rate, split="\\s+")[[1]][3])
        hours.remaining <- mAh / mA
    }
    else stop("Invalid 'state' string:", state)
    list(time=date(), hours.remaining=hours.remaining)
}

count.lines <- function(file) wc()["lines"]

memory.usage <- function(n, what=c("double","integer")) {
    what <- match.arg(what)
    nbytes <- c(integer=4, double=8)
    c(Gb=n * nbytes[what] * 2^(-30))
}

extract.lines <- function(lines, fromfile, tofile, method=c("sed","line")) {
    switch(match.arg(method),
           "sed" = {
               call <- paste("-e", paste(lines, "p", sep="", collapse=" -e "))
               call <- paste("sed -n", call, fromfile, ">", tofile)
           },
           "line" = {
               call <- paste("line", paste(lines, collapse=","), fromfile, ">", tofile)
           })
    print(call)
    system(call)
}

file.strsub <- function(str, sub, dir=".") {
    ## Replace string 'str' with 'sub' in all filenames in 'dir'

    f <- dir(dir)
    ff <- gsub(str, sub, f)
    for(i in seq(along=f))
        file.rename(paste(dir, f[i], sep="/"), paste(dir, ff[i], sep="/"))
}
## block-7 ends here

## [[file:~/src/dedR/ded.R.org::*Graphics][block-8]]
plot.densCols <- function(x, y=NULL, palette="Blues", ...) {
    require(geneplotter)
    require(RColorBrewer)
    stopifnot(palette %in%
              c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys",
                "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples",
                "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd",
                "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy",
                "RdYlBu", "RdYlGn", "Spectral"))
    plot(x, y, col=densCols(x, y, colramp=colorRampPalette(brewer.pal(9, palette)[-(1:3)])), ...)
}

plot.blank <- function(xlim=c(0,1), ylim=c(0,1), xlab="", ylab="")
    plot(NA, xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", bty="n", type="n", xlab=xlab, ylab=ylab)

my.brewer.pal <- function(n, type=c("qualitative","sequential","divergent"), name) {
    require(RColorBrewer)
    switch(match.arg(type),
           qualitative = {
               stopifnot(n <= 9)
               idx <- (1:9) ## don't
               if(n < 9) idx <- idx[-6][1:n]
               if(missing(name)) name <- "Set1"
               brewer.pal(9, name)[idx]
           },
           sequential = {
               if(missing(n)) n <- 9
               if(missing(name)) name <- "RdPu" ## "BuGn" "Reds" etc
               brewer.pal(n, name)
           },
           divergent = stop("unimplemented"))
}

zoomplot <- function(x, y, ...) {
    plot(x, y, ...)
    while(TRUE) {
        newx <- locator(2)$x
        if(newx[1] == newx[2]) break
        plot(x, y, xlim=newx, type="n", ...) ## plot invisibly to establish (nice) axis limits
        axp <- par("xaxp")
        axgrid <- seq(axp[1], axp[2], by=0.1) ## obtain desired axis ticks
        axp <- c(range(axgrid), length(axgrid)-1)
        plot(x, y, xlim=newx, xaxp=axp, ...) ## and plot visibly with desired axis ticks
        abline(v=axgrid, col="lightgrey") ## draw vertical lines there
    }
}


gradient.under.graph <- function(n=100, y0=0, mu=-0.7, sd=1, nrects=1000) {
    y <- y0 + cumsum(rnorm(n, mean=mu, sd=sd))
    plot(NA, xlim=c(1,n), ylim=range(y), bty="n")
    col <- colorRampPalette(colors=c("dark blue","light blue"))(nrects)
    incr <- (max(y) - min(y)) / nrects
    rect(0, seq(min(y), max(y)-incr, length=nrects), n, seq(min(y)+incr, max(y), length=nrects), col=col, border=NA)
    polygon(x=c(1,1:n,n), y=c(max(y), y, max(y)), col="white", border=NULL)
}

print.binary.matrix <- function(X) {
    dimn <- dim(X)
    stopifnot(length(dimn) %in% c(2,3))
    dimn3 <- if(length(dimn) == 3) dimn[3] else 1

    for(k in 1:dimn3) {
        if(dimn3 > 1) cat(dimnames(X)[[3]][k], "\n") ## OK even if dimnames are NULL
        for(i in 1:(dimn[1])) {
            for(j in 1:(dimn[2])) {
                if(dimn3 > 1) {
                    if(is.na(X[i,j,k])) cat(" ")
                    else if(X[i,j,k]) cat(if(class(X)[1] == "haplotypes") "1" else "1 ")
                    else cat(if(class(X)[1] == "haplotypes") "." else ". ")
                }
                else {
                    if(is.na(X[i,j])) cat(" ")
                    else if(X[i,j]) cat("1")
                    else cat(".")
                }
            }
            cat("\n")
        }
        cat("\n")
    }
}
## block-8 ends here

## [[file:~/src/dedR/ded.R.org::*Graphics][block-9]]
## plot a bunch of KDEs contained in a list
plot.list.of.densities <- function(dd, col, lty, lwd, xlab="", ylab="", xmin, xmax, legend=FALSE,
                                   add=FALSE, ...) {
    ymax <- max(sapply(dd, function(d) max(d$y)), na.rm=TRUE)
    if(missing(xmax)) xmax <- max(sapply(dd, function(d) max(d$x)), na.rm=TRUE)
    if(missing(xmin)) xmin <- min(sapply(dd, function(d) min(d$x)), na.rm=TRUE)

    if(!add) plot(x=c(xmin, xmax), y=c(0, ymax), type="n", xlab=xlab, ylab=ylab, ...)
    n <- length(dd)
    if(missing(col)) col <- 1:n else col <- rep(col, length=n)
    if(missing(lty)) lty <- rep(1, n) else lty <- rep(lty, length=n)
    if(missing(lwd)) lwd <- rep(1, n) else lwd <- rep(lwd, length=n)

    if(legend) {
        stopifnot(!is.null(names(dd)))
        legend("topright", legend=names(dd), col=col, lty=lty, lwd=lwd, bty="n")
    }
    for(i in 1:n)
        lines(dd[[i]], col=col[i], lty=lty[i], lwd=lwd[i])
}
## block-9 ends here

## [[file:~/src/dedR/ded.R.org::*Plot%20columns][block-10]]
plot.columns <- function(A, B=NULL, densCols=FALSE, lm=FALSE, print.corr=FALSE, label=TRUE,
                         names=NULL, main=NULL, col=1, adj=10, ...) {
    if(is.null(B)) {
        symmetric = TRUE
        B <- A
    }
    else symmetric=FALSE
    stopifnot(nrow(A) == nrow(B))
    if(densCols) require("prada")
    nA <- ncol(A) ; nB <- ncol(B)
    op <- par(mfrow=c(nA, nB), mar=rep(0,4), oma=if(is.null(main)) c(0,0,0,0) else c(0,0,3,0)) ; on.exit(par(op))
    if(label && is.null(names)) names <- c("A","B")
    names <- rev(names)
    for(i in 1:nA) {
        for(j in 1:nB) {
            x <- B[,j]
            y <- A[,i]
            if(symmetric && j == i) {
                plot.blank()
                if(!is.null(colnames(A)))
                    text(mean(range(x)), mean(range(y)), colnames(A)[i], cex=1.2)
                next
            }
            if(densCols) col <- densCols(x,y)
            plot(x, y, col=col, xaxt="n", yaxt="n", mgp=c(.1,.1,.1), ...)
            if(lm) abline(lm(y ~ x), col="red")
            if(print.corr) {
                r <- round(cor(y, x), 10)
                legend("bottomleft", legend=r, bty="n", bg=NA, text.col="red", cex=.5+r)
            }
            if(label) {
                xmid <- mean(range(x, na.rm=TRUE))
                xmin <- min(x, na.rm=TRUE)
                xmin <- xmin + 2*(xmid - xmin) / adj
                ymid <- mean(range(y, na.rm=TRUE))
                ymin <- min(y, na.rm=TRUE)
                ymin <- ymin + (ymid - ymin) / adj
                text(x=c(xmid, xmin), y=c(ymin, ymid), labels=paste(names, c(j,i), sep=""), col="red", cex=1.5)
            }
            ## text(x=1.1 * min(range(A[,i], na.rm=TRUE)), y=1.1 * min(range(B[,j], na.rm=TRUE)),
            ## labels=round(cor(A[,i], B[,j]), 2), col="red", cex=1.5)
            ## plot.blank()
            ## cell <- c(i,j)
            ## text(x=0.5, y=0.5, paste(cell[1], cell[2]), cex=2)
        }
    }
    if(!is.null(main)) mtext(text=main, outer=TRUE, cex.main=list(...)$cex.main)
}
## block-10 ends here

## [[file:~/src/dedR/ded.R.org::*Colours][block-11]]
rgb.range <- function(x, col=c("red","green","blue")) {
    ok <- !is.na(x)
    x <- x / max(x, na.rm=TRUE)

    rgb <- c("red","green","blue")
    cols <- replicate(3, rep(1, length(x)), simplify=FALSE)
    names(cols) <- rgb
    col <- match.arg(col)
    for(c in rgb) {
        if(c == col) cols[[c]][ok] <- x[ok]
        else cols[[c]][ok] <- 0
    }
    do.call("rgb", cols)
}

## to convert output of color ramp functions
rgb.matrix <- function(x) {
    dimnames(x) <- list(NULL, c("red","green","blue"))
    do.call("rgb", c(as.list(as.data.frame(x)), maxColorValue=255))
}

my.pars <- function() list(bg.black=list(bg="black", col.axis="lightgrey", col.lab="lightgrey",
                           col.main="lightgrey", fg="lightgrey"),
                           vanilla=list(xaxt="n", yaxt="n", bty="n", mar=rep(0,4)))
## block-11 ends here

## [[file:~/src/dedR/ded.R.org::*Stats][block-12]]
rdirichlet <- function(n, a) {
    ## n draws from Dirichlet(a_1, ..., a_K)
    K <- length(a)
    x <- matrix(rgamma(n*K, shape=a, rate=1), nrow=n, ncol=K, byrow=TRUE)
    x / rowSums(x)
}


apriori <- function(x, t) {
    ## "Apriori" algorithm to form groups of variables (columns of
    ## x). Based on description in EoSL (Hastie, Tibshirani et al.)

    stop("Not correct or even vaguely properly thought-out")
    
    apriori.print.clusters <- function() {
        for(g in clusters)
            cat(paste(g, "\n"))
    }
    
    singletons <- which(colMeans(x) > t)
    k <- 1
    clusters <- as.vector(singletons, mode="list")
    while(TRUE) {
        ## clusters <- newclusters
        ## newclusters <- NULL
        k <- k+1
        cat("k =", k, "\n")
        for(i in singletons) {
            cat(i, "\n")
            goti <- FALSE
            recipients <- which(sapply(clusters, length) == k-1)
            if(k == 2) recipients <- recipients[recipients != i]
            for(j in recipients) {
                g <- sort(c(clusters[[j]], i))
                stopifnot(length(g) == k)
                if(any(sapply(clusters, function(gg) identical(gg, g))))
                    next
                support <- mean(rowSums(x[,g,drop=FALSE]) == k)
                if(support > t) {
                    clusters <- c(clusters[-j], list(g))
                    goti <- TRUE
                }
            }
            apriori.print.clusters()
            if(goti) singletons <- singletons[singletons != i]
        }
    }
    clusters
}

kde2D <- function(x, bandwidth, gridsize) {
    ## default bandwidth taken from .smoothScatterCalcDensity
    if (missing(bandwidth))
        bandwidth <- diff(apply(x, 2,
                                quantile, probs = c(0.05, 0.95), na.rm = TRUE))/25
    if (missing(gridsize)) gridsize <- c(51,51)
    else if (length(gridsize) == 1)
        gridsize <- rep(gridsize, 2)
    
    kde <- bkde2D(x, bandwidth=bandwidth, gridsize=gridsize)
    structure(kde$fhat, dimnames=list(kde$x1, kde$x2))
}


pvalue.lm <- function(lm.fit) {
    f <- summary.lm(lm.fit)$fstatistic
    if(is.null(f)) f <- rep(NA, 3)
    f <- as.list(f)
    names(f) <- c("q","df1","df2")
    do.call("pf", c(f, lower.tail=FALSE))
}

anova.pvalue.table <- function(anovas) {
    ## Construct a table of p-values from a list of ANOVA tables
    rnn <- sapply(anovas, rownames)
    stopifnot(apply(rnn, 1, all.same))
    ptab <- sapply(anovas, "[[", "p-value")
    rownames(ptab) <- rnn[,1]
    ptab
}

estimate.mode <- function(sample, discrete=FALSE) {
    if(discrete) as.integer(names(which.max(table(sample)))) ## dodgy?
    else {
        d <- density(sample, n=1024)
        d$x[which.max(d$y)]
    }
}

mcmc.plot <- function(x, from, to, prior=NULL) {
    ## I think coda::plot.mcmc should have the axes reversed so things
    ## line up
    if(is.null(dim(x))) dim(x) <- c(length(x),1)
    p <- ncol(x)
    par(mfcol=(c(2,p)))
    for(i in 1:p) {
        ## plot(density(x[,i], from=from, to=to), bty="n", yaxt="n", xlab="", ylab="", main="")
        hist(x[,i], bty="n", yaxt="n", xlab="", ylab="", main="")
        if(!is.null(prior)) curve(prior[[i]](x), col="grey", add=TRUE)
        plot(x[,i], seq(nrow(x)), xlim=c(from,to), type="l", bty="l", xlab=colnames(x)[i], ylab="")
    }
}


qqgen.multi <- function(pp, f=identity, drawfun=list(points), add=FALSE, col="black", lwd=1,...){
    ylim <- range(f(unlist(pp)), na.rm=TRUE)
    nmax <- max(sapply(pp, length))
    xlim <- range(f((1:nmax)/(nmax+1)))

    col <- rep(col, length.out=length(pp))
    lwd <- rep(lwd, length.out=length(pp))
    drawfun <- rep(drawfun, length.out=length(pp))

    for(i in 1:length(pp))
        qqgen(pp[[i]], f=f, drawfun=drawfun[[i]], add=(i>1 || add), col=col[i], lwd=lwd[i], xlim=xlim, ylim=ylim, ...)
}

qqgen <- function(p, f=identity, add=FALSE, drawfun=points, plotci=FALSE, ...) {
    ## use f to transform to -log10 scale, or to test statistics via q{dist}
    o <- f(sort(p)) ## removes NAs
    o <- o[is.finite(o)]
    n <- length(o)
    e <- f((1:n)/(n+1))
    if(!add) {
        plot(e, o, type="n", ...)
        abline(c(0,1), col="black")
    }
    drawfun(e, o, ...)
}


qq.plotci <- function(e, f, df=2) {
    stop("Need to work out how to use this fragment taken from plot.qqchisq.for.paper")
    n <- length(e)
    out <- seq(1, n, length.out=10000)
    back <- seq(n, 1, length.out=10000)
    upper <- qchisq(qbeta(0.975, back, out), df=df, lower.tail=FALSE)
    lower <- qchisq(qbeta(0.025, back, out), df=df, lower.tail=FALSE)
    polygon(c(e[out], e[back]), c(lower, rev(upper)), density=NA, col="grey")
}

qqchisq <- function(pvals, df, new=TRUE, lim=NULL, cols=NULL, simulate.from.null=FALSE, ndraws,
                    hack=FALSE, ...) {
    if(any(is.na(pvals))) stop("There are NAs: is this code correct in that case?")
    n <- length(pvals)
    ord <- order(1-pvals)
    if(is.null(cols)) cols <- rep("black", n)
    else if(length(cols) == 1) cols <- rep(cols, n)
    else if(is.numeric(cols) && length(cols) == n) cols <- rgb.range(cols, col="red")
    else if(length(cols) != n) stop("Having problems interpreting colour argument")
    obs.q <- qchisq(1 - pvals, df=df)
    exp.q <- qchisq(1:n/(n+1), df=df)
    if(new && is.null(lim)) {
        vals <- c(exp.q, obs.q)
        vals[is.infinite(vals)] <- NA
        lim <- range(vals, na.rm=TRUE)
    }
    qq <- qqplot(exp.q, obs.q, type="n", xlim=lim, ylim=lim, plot.it=new, ...)

    if(FALSE && hack) {
        red <- cols == "red"
        do.call("points", list(x=qq$x[red], y=qq$y[red], col="red", ...))
    }

    if(simulate.from.null) {
        robs.q <- replicate(ndraws, rchisq(n, df=df), simplify=FALSE)
        rqq <- lapply(robs.q, function(q) qqplot(exp.q, q, plot.it=FALSE))
        sapply(rqq, function(q) lines(q$x, q$y, col="grey"))
    }
    abline(0, 1, col="grey" )
    do.call("points", list(x=qq$x, y=qq$y, col=cols[ord], ...))
    invisible(qq)
}

rqqchisq <- function(pvals, df, cols=NULL, simulate.from.null=FALSE, ndraws, ...) {
    qq <- replicate(ndraws, qqchisq(1:N, df, simulate.from.null=TRUE), simplify=FALSE)
    cols <- rainbow(ndraws)
    plot(qq[[1]]$x, qq[[1]]$y, col=cols[1], type="l", ...)
    sapply(2:ndraws, function(draw) lines(qq[[draw]]$x, qq[[draw]]$y, col=cols[draw]))
    abline(c=c(0,1))
    invisible(qq)
}
## block-12 ends here

## [[file:~/src/dedR/ded.R.org::*Partial%20eigen][block-13]]
## http://www.netlib.org/lapack/explore-html/dsyevr.f.html
## http://tolstoy.newcastle.edu.au/R/devel/04/11/1258.html
## (R-devel posts from Jon McAuliffe and Brian Ripley)
partial.eigen <- function(x, only.values=FALSE, nvectors=ncol(x), tol=1e-5, LWORK, LIWORK=10*N) {
    ## Compute first 'nvectors' eigenvectors and their eigenvalues of x.
    ## The fortran function is supposed to compute the optimal sizes LWORK/LIWORK of the double/int
    ## 'work areas' WORK/IWORK, by supplying -1 as LWORK/LIWORK. However, I can only get this to work
    ## for the work space of doubles LWORK.
    ## 2007.02.10

    ## I think that's because you didn't do it carefully enough; both are returned in the first element of
    ## respective arrays; one array is double so you have to cast that to int

    stopifnot(nrow(x) == ncol(x))
    N <- ncol(x)

    args <- list(JOBZ = if(only.values) 'N' else 'V',
                 RANGE = 'I',                                                        # specify desired eigenvectors with IL, IU args
                 UPLO = 'L',                                                         # it uses the values in the lower triangle
                 N = as.integer(N),                                                  # the 'order' of the matrix
                 A = as.double(x),
                 LDA = as.integer(N),                                                # increment to get from one column to next
                 VL = as.double(-1),                                                 # use these if you want all eigenvectors whose
                 VU = as.double(-1),                                                 # evalues fall within (VL, VU]
                 IL = as.integer(N - nvectors + 1),
                 IU = as.integer(N),
                 ABSTOL=as.double(tol),
                 M = integer(1),
                 W = double(nvectors),                                               # space to store evalues
                 Z = double(N * nvectors),                                           # space to store evectors
                 LDZ = as.integer(N),                                                #
                 ISUPPZ = integer(2*nvectors),

                 WORK = NA,                                                          # 'work area' of doubles
                 LWORK = NA,                                                         # size of 'work area' of doubles

                 ## IWORK = if(compute.LIWORK) double(1) else double(LIWORK),        # 'work area' of ints
                 ## 29 Nov 2007 why are you passing double rather than int here?
                 IWORK = NA,           # 'work area' of ints
                 LIWORK = NA, # size of 'work area' of ints

                 INFO=integer(1))

    ## First request optimal area of work arrays
    args$WORK <- double(1)
    args$LWORK <- as.integer(-1)
    args$IWORK <- integer(1)
    args$LIWORK <- as.integer(-1)

    out <- do.call(".Fortran", c(list("dsyevr"), args))

    cat("LWORK =", out$WORK, "\nLIWORK =", out$IWORK, "\n")

    ## Now do real thing
    args$WORK <- double(args$LWORK)
    args$LWORK <- as.integer(args$LWORK)
    args$IWORK <- integer(args$LIWORK)
    args$LIWORK <- as.integer(args$LIWORK)

    out <- do.call(".Fortran", c(list("dsyevr"), args))

    list(values=rev(out$W[1:nvectors]), vectors=matrix(out$Z,ncol=nvectors)[,nvectors:1])
}
## block-13 ends here

## [[file:~/src/dedR/ded.R.org::*Partial%20eigen%20orig][block-14]]
## http://www.netlib.org/lapack/explore-html/dsyevr.f.html
## http://tolstoy.newcastle.edu.au/R/devel/04/11/1258.html
## (R-devel posts from Jon McAuliffe and Brian Ripley)
partial.eigen.orig <- function(x, only.values=FALSE, nvectors=ncol(x), tol=1e-5,
                               compute.LWORK=FALSE,
                               LWORK, compute.LIWORK=FALSE, LIWORK=10*N) {
    ## Compute first 'nvectors' eigenvectors and their eigenvalues of x.
    ## The fortran function is supposed to compute the optimal sizes LWORK/LIWORK of the double/int
    ## 'work areas' WORK/IWORK, by supplying -1 as LWORK/LIWORK. However, I can only get this to work
    ## for the work space of doubles LWORK.
    ## 2007.02.10

    ## I think that's because you didn't do it carefully enough; both are returned in the first element of
    ## respective arrays; one array is double so you have to cast that to int




    stopifnot(nrow(x) == ncol(x))
    N <- ncol(x)

    out <- .Fortran("dsyevr",
                    JOBZ = if(only.values) 'N' else 'V',
                    RANGE = 'I',                                                        # specify desired eigenvectors with IL, IU args
                    UPLO = 'L',                                                         # it uses the values in the lower triangle
                    N = as.integer(N),                                                  # the 'order' of the matrix
                    A = as.double(x),
                    LDA = as.integer(N),                                                # increment to get from one column to next
                    VL = as.double(-1),                                                 # use these if you want all eigenvectors whose
                    VU = as.double(-1),                                                 # evalues fall within (VL, VU]
                    IL = as.integer(N - nvectors + 1),
                    IU = as.integer(N),
                    ABSTOL=as.double(tol),
                    M = integer(1),
                    W = double(nvectors),                                               # space to store evalues
                    Z = double(N * nvectors),                                           # space to store evectors
                    LDZ = as.integer(N),                                                #
                    ISUPPZ = integer(2*nvectors),
                    WORK = if(compute.LWORK) double(1) else double(LWORK),              # 'work area' of doubles
                    LWORK = if(compute.LWORK) as.integer(-1) else as.integer(LWORK),    # size of 'work area' of doubles

                    ## IWORK = if(compute.LIWORK) double(1) else double(LIWORK),           # 'work area' of ints
                    ## 29 Nov 2007 why are you passing double rather than int here?
                    IWORK = if(compute.LIWORK) integer(1) else integer(LIWORK),           # 'work area' of ints

                    LIWORK = if(compute.LIWORK) as.integer(-1) else as.integer(LIWORK), # size of 'work area' of ints

                    INFO=integer(1))

    if(compute.LWORK || compute.LIWORK) lapply(out[c("WORK","IWORK")], "[", 1)
    else list(values=rev(out$W[1:nvectors]), vectors=matrix(out$Z,ncol=nvectors)[,nvectors:1])
}
## block-14 ends here

## [[file:~/src/dedR/ded.R.org::*My%20eigen][block-15]]
my.eigen <- function (x, symmetric, only.values = FALSE, EISPACK = FALSE)
{
    x <- as.matrix(x)
    if (!is.null(dimnames(x)))
        dimnames(x) <- list(NULL, NULL)
    n <- nrow(x)
    if (!n)
        stop("0 x 0 matrix")
    if (n != ncol(x))
        stop("non-square matrix in 'eigen'")

    complex.x <- is.complex(x)
    if (!complex.x && !is.numeric(x))
        stop("numeric or complex values required in 'eigen'")

    if (missing(symmetric))
        symmetric <- isSymmetric.matrix(x)
    will.use <- if(symmetric) lower.tri(x, diag=TRUE) else TRUE
    if (any(!is.finite(x[will.use])))
        stop("infinite or missing values in 'x'")
    if (is.numeric(x))
        storage.mode(x) <- "double"
    if (!EISPACK) {
        if (symmetric) {
            z <- if (!complex.x)
                .Call("La_rs", x, only.values, PACKAGE = "base")
            else .Call("La_rs_cmplx", x, only.values, PACKAGE = "base")
            ord <- rev(seq_along(z$values))
        }
        else {
            z <- if (!complex.x)
                .Call("La_rg", x, only.values, PACKAGE = "base")
            else .Call("La_rg_cmplx", x, only.values, PACKAGE = "base")
            ord <- sort.list(Mod(z$values), decreasing = TRUE)
        }
        return(list(values = z$values[ord], vectors = if (!only.values) z$vectors[,
                                            ord, drop = FALSE]))
    }
    dbl.n <- double(n)
    if (symmetric) {
        if (complex.x) {
            xr <- Re(x)
            xi <- Im(x)
            z <- .Fortran("ch", n, n, xr, xi, values = dbl.n,
                          !only.values, vectors = xr, ivectors = xi, dbl.n,
                          dbl.n, double(2 * n), ierr = integer(1), PACKAGE = "base")
            if (z$ierr)
                stop(gettextf("'ch' returned code %d in 'eigen'",
                              z$ierr), domain = NA)
            if (!only.values)
                z$vectors <- matrix(complex(real = z$vectors,
                                            imaginary = z$ivectors), ncol = n)
        }
        else {
            z <- .Fortran("rs", n, n, x, values = dbl.n, !only.values,
                          vectors = x, dbl.n, dbl.n, ierr = integer(1),
                          PACKAGE = "base")
            if (z$ierr)
                stop(gettextf("'rs' returned code %d in 'eigen'",
                              z$ierr), domain = NA)
        }
        ord <- sort.list(z$values, decreasing = TRUE)
    }
    else {
        if (complex.x) {
            xr <- Re(x)
            xi <- Im(x)
            z <- .Fortran("cg", n, n, xr, xi, values = dbl.n,
                          ivalues = dbl.n, !only.values, vectors = xr,
                          ivectors = xi, dbl.n, dbl.n, dbl.n, ierr = integer(1),
                          PACKAGE = "base")
            if (z$ierr)
                stop(gettextf("'cg' returned code %d in 'eigen'",
                              z$ierr), domain = NA)
            z$values <- complex(real = z$values, imaginary = z$ivalues)
            if (!only.values)
                z$vectors <- matrix(complex(real = z$vectors,
                                            imaginary = z$ivectors), ncol = n)
        }
        else {
            z <- .Fortran("rg", n, n, x, values = dbl.n, ivalues = dbl.n,
                          !only.values, vectors = x, integer(n), dbl.n,
                          ierr = integer(1), PACKAGE = "base")
            if (z$ierr)
                stop(gettextf("'rg' returned code %d in 'eigen'",
                              z$ierr), domain = NA)
            ind <- z$ivalues > 0
            if (any(ind)) {
                ind <- seq.int(n)[ind]
                z$values <- complex(real = z$values, imaginary = z$ivalues)
                if (!only.values) {
                    z$vectors[, ind] <- complex(real = z$vectors[,
                                                ind], imaginary = z$vectors[, ind + 1])
                    z$vectors[, ind + 1] <- Conj(z$vectors[, ind])
                }
            }
        }
        ord <- sort.list(Mod(z$values), decreasing = TRUE)
    }
    list(values = z$values[ord], vectors = if (!only.values) z$vectors[,
                                 ord, drop = FALSE])
}
## block-15 ends here

## [[file:~/src/dedR/ded.R.org::*Org%20mode][block-16]]
source.babel <- function(file)
    source(pipe(paste("~/bin/org-babel-tangle R", file)))

array2org <- function(x) {

    if(is.null(cn <- colnames(x)))
        cn <- seq(ncol(x))

    cat("|", paste(cn, collapse=" | "), "|\n")
    cat(paste("|", paste(rep("--", ncol(x)), collapse="+"), "|", sep=""), "\n")
    for(i in 1:nrow(x))
        cat("|", paste(x[i,], collapse=" | "), "|\n")
}
## block-16 ends here

## [[file:~/src/dedR/ded.R.org::*Genomics][block-17]]
query.ucsc.rs.numbers <- function(rs, database="hg18") {
    table <- c(hg17="snp125", hg18="snp130")
    condition <- paste("name=", "'", rs, "'", sep="", collapse=" OR ")
    query.ucsc(fields=c("name", "chrom", "chromStart", "chromEnd"), condition=condition, database=database, table=table[database])
}

query.ucsc <- function(fields, condition, database="hg18", table="snp130") {
    fields <- paste(fields, collapse=",")
    call <- "mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A"
    call <- paste(call, "--database=", database, "-B -e ")
    call <- paste(call, '"', 'SELECT ', fields, 'FROM', table, 'WHERE ', '(', condition, ')', '"')
    cat(call, "\n")
    read.delim(pipe(call), as.is=TRUE)
}
## block-17 ends here
