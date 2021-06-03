# This code is taken from rpart.plot package
# Originally written by Stephen Milborrow (2010) and
# derived from the rpart code by written by
#   Terry M Therneau and Beth Atkinson:
#   http://mayoresearch.mayo.edu/mayo/research/biostat/splusfunctions.cfm
# and the R port and modifications of that code by Brian Ripley:
#   www.stats.ox.ac.uk/~ripley/
#
# The functions were modified by Anna Neufeld in 9/2020.
# The only changes from the rpart.plot package are that pvalues and confidence intervals
# are now printed.
#------------------------------------------------------------------------------

# allowable values of prp's type argument
TYPE0.default         <- 0
TYPE1.all             <- 1
TYPE2.all.under       <- 2
TYPE3.fancy.no.all    <- 3
TYPE4.fancy.all       <- 4
TYPE5.varname.in.node <- 5


#' A copy of the workhorse function from rpart.plot that has been modified.
#' @keywords internal
#' @noRd
inner.plot <- function(x=stop("no 'x' arg"),
    type=2, extra="auto", under=FALSE, fallen.leaves=TRUE,
    digits=2, varlen=0, faclen=0, roundint=TRUE,
    cex=NULL, tweak=1,
    clip.facs=FALSE, clip.right.labs=TRUE,
    snip=FALSE,
    box.palette="auto", shadow.col=0,
    ...)
{
 if(!inherits(x, "rpart"))
        stop("Not an rpart object")
    # We have to get "trace" for get.modelframe.info, but I don't want to
    # add trace to the rpart.plot arg list, hence the following bit of
    # code to get trace from the dots.
    dots <- match.call(expand.dots=FALSE)$...
    trace <- 0
    if(!is.null(dots$trace))
        trace <- eval(dots$trace)
    x$varinfo <- rpart.plot:::get.modelframe.info(x, roundint, trace,
                                     parent.frame(), "rpart.plot")
    prp(x,
        type=type, extra=extra, under=under, fallen.leaves=fallen.leaves,
        digits=digits, varlen=varlen, faclen=faclen, roundint=roundint,
        cex=cex, tweak=tweak,
        clip.facs=clip.facs, clip.right.labs=clip.right.labs,
        snip=snip,
        box.palette=box.palette, shadow.col=shadow.col,
        ...)
}

#' #' @keywords internal
#' #' @noRd
#' rpart.plot.version1 <- function(x=stop("no 'x' arg"),
#'     type=0, extra=0, under=FALSE, fallen.leaves=FALSE,
#'     digits=2, varlen=-8, faclen=3,
#'     cex=NULL, tweak=1,
#'     snip=FALSE,
#'     box.palette=0, shadow.col=0,
#'     ...)
#' {
#'     if(!inherits(x, "rpart"))
#'         stop("Not an rpart object")
#'     # We have to get "trace" for get.modelframe.info, but I don't want to
#'     # add trace to the rpart.plot arg list, hence the following bit of
#'     # code to get trace from the dots.
#'     dots <- match.call(expand.dots=FALSE)$...
#'     trace <- 0
#'     if(!is.null(dots$trace))
#'         trace <- eval(dots$trace)
#'     ## prepared to fixin 1 min
#'     x$varinfo <- rpart.plot:::get.modelframe.info(x, roundint=FALSE, trace,
#'                                      parent.frame(), "rpart.plot.version1")
#'     prp(x,
#'         type=type, extra=extra,
#'         under=under, fallen.leaves=fallen.leaves, clip.facs=FALSE,
#'         digits=digits, varlen=varlen, faclen=faclen, roundint=FALSE,
#'         cex=cex, tweak=tweak,
#'         snip=snip,
#'         box.palette=box.palette, shadow.col=shadow.col,
#'         ...)
#' }

#' @keywords internal
#' @noRd
#' @importFrom graphics axis grid par plot rect text
prp <- function(x=stop("no 'x' arg"),
    type=0, extra=0, under=FALSE, fallen.leaves=FALSE,
    nn=FALSE, ni=FALSE, yesno=TRUE,
    branch=if(fallen.leaves) 1 else .2,
    uniform=TRUE, left=TRUE, xflip=FALSE, yflip=FALSE,
    digits=2, varlen=-8, faclen=3, roundint=TRUE,
    cex=NULL, tweak=1,
    clip.facs=FALSE, clip.right.labs=TRUE,
    compress=TRUE, ycompress=uniform,
    Margin=0, space=1, gap=NULL,
    snip=FALSE, snip.fun=NULL, trace=FALSE,

    box.col=0, box.palette=0,
    pal.thresh=NULL, pal.node.fun=FALSE,
    border.col=col,
    round=NULL, leaf.round=NULL,
    shadow.col=0, prefix="", suffix="", xsep=NULL,

    under.percent=2, under.font=font, under.col=1, under.cex=.8,

    split.cex=1, split.font=2, split.family=family, split.col=1,
    split.box.col=0, split.border.col=0,
    split.lty=1, split.lwd=NULL, split.round=0,
    split.shadow.col=0,
    split.prefix="", right.split.prefix=NULL,
    split.suffix="", right.split.suffix=NULL,
    facsep=",", eq=" = ", lt=" < ", ge=" >= ",

    branch.col=if(rpart.plot:::is.zero(branch.type)) 1 else "gray",
    branch.lty=1, branch.lwd=NULL,
    branch.type=0, branch.tweak=1,
    min.branch.width=.002, branch.fill=branch.col,

    nn.cex=NULL, nn.font=3, nn.family="", nn.col=1,
    nn.box.col=0, nn.border.col=nn.col,
    nn.lty=1, nn.lwd=NULL, nn.round=.3,
    yes.text="yes", no.text="no",

    node.fun=NULL,
    split.fun=NULL,
    FUN="text",

    nspace=branch, minbranch=.3, do.par=TRUE,
    add.labs=TRUE,
    clip.left.labs=(type == 5),
    fam.main="",
    yshift=0, yspace=space, shadow.offset=.4,

    split.adj=NULL, split.yshift=0, split.space=space,
    split.yspace=yspace, split.shadow.offset=shadow.offset,

    nn.adj=.5, nn.yshift=0, nn.space=.8, nn.yspace=.5,

    ygap=gap/2, under.ygap=.5, yesno.yshift=0,
    xcompact=TRUE, ycompact=uniform, xcompact.ratio=.8, min.inter.height=4,

    max.auto.cex=1, min.auto.cex=.15, ycompress.cex=.7, accept.cex=1.1,
    shift.amounts=c(1.5, 2),
    Fallen.yspace=.1,
    boxes.include.gap=FALSE,
    legend.x=NULL, legend.y=NULL, legend.cex=1,
    ...)
{
    check.dots <- function(dots) # check dots arguments, if any
    {
        legal.dots.args <- # they are legal if we have code to process them later
            c("adj", "cex.main", "cex.sub", "col", "col.main", "col.sub", "family",
              "font", "lty", "lwd", "main", "mar", "sub", "xlim", "xpd", "ylim")

        if(length(dots) > 0) {
            names <- names(dots)
            pmatch <- pmatch(names, legal.dots.args, duplicates.ok=TRUE)
            if(any(is.na(pmatch))) {
                # report the first illegal arg
                ibad <- (1:length(dots))[is.na(pmatch)]
                stop("prp: illegal argument \"", names[ibad][1], "\"")
            }
            duplicated <- duplicated(pmatch)
            if(any(duplicated))
                stop("prp: duplicated argument \"", names[duplicated][1], "\"")
        }
    }
    merge1 <- function(vec, split.vec)
    {
        split.vec <-rpart.plot:::recycle(split.vec, nodes)
        split.vec[is.leaf] <-rpart.plot:::recycle(vec, nodes)[is.leaf]
        split.vec
    }
    draw.labs <- function(draw.shadows1, draw.split.shadows1)
    {
        # put the labels on the screen, text after \n\n (if any) goes under the box
        draw.labs1 <- function(labs, boxes, yspace, cex, font, family, col,
                               draw.shadows1, make.space.for.shadows, shadow.col, round)
        {
            draw.under.text <- function() { # draw the text under the box, and its shadow
                height1 <- rpart.plot:::my.strheight("M", cex, font, family)
                cex <- under.cex * cex
                under.height <- rpart.plot:::my.strheight(sep.labs$under.box, cex, under.font, family)
                x <- xy$x
                y <- boxes$y1 - under.ygap * height1 - .5 * under.height
                width  <- .5 * rpart.plot:::my.strwidth(sep.labs$under.box,  cex, under.font, family)
                height <- .5 * rpart.plot:::my.strheight(sep.labs$under.box, cex, under.font, family)
                # the magic numbers 1.4 and 1.2 seem about right visually
                if(make.space.for.shadows)
                    height <- 1.4 * height
                if(draw.shadows1)
                    rpart.plot:::draw.shadow(x - 1.2 * width, y - height,
                         x + 1.2 * width, y + height,
                         xlim, ylim, 0, shadow.col, shadow.offset)
                else {
                    # draw a white box to hide anything underneath, then draw the text
                    rect(x - 1.2 * width, y - height,
                         x + 1.2 * width, y + height, col=bg, border=0)
                    FUN(x, y, sep.labs$under.box,
                        cex=cex, font=under.font, family=family, col=under.col)
                }
            }
            #--- draw.labs1 starts here ---
            FUN <- rpart.plot:::check.func.args(FUN, "FUN argument to the prp", graphics::text)
            sep.labs <- separate.labs(labs)
            xy <- get.box.centers(boxes)
            # draw the text after \n\n if any under the box
            if(!all(nchar(sep.labs$under.box) == 0))
                draw.under.text()
            # draw the text before \n\n in the box
            if(!draw.shadows1)
                FUN(xy$x, xy$y, sep.labs$in.box,
                    cex=cex, font=font, family=family, col=col)
        }
        #--- draw.labs starts here ---
        if(boxes.include.gap) {
            # For debugging: make get.boxes expand the boxes to include
            # what would normally be the gap between the boxes.
            # With optimum cex, at least one pair of boxes will just touch.
            print("boxes.include.gap is TRUE\n")
            split.space <- split.space + gap/2
            split.yspace <- split.yspace + ygap/2
            space <- space + gap/2
            yspace <- yspace + ygap/2
            gap <- ygap <- 0
        }
        # With type==TYPE2.all.under and no visible split box, branch lines
        # look better if just a small space around labs.
        small.underspace <- type == TYPE2.all.under &&
            rpart.plot:::is.box.invisible(split.box.col, split.border.col, bg)

        split.boxes <-
            draw.boxes(if(is.fancy(type)) "left" else "default", draw.split.shadows1,
                   split.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   split.cex * cex, split.font, split.family, split.adj, split.yshift,
                   split.box.col, split.border.col,
                   split.lty, split.lwd, split.space, split.yspace,
                   split.round * split.strheight,
                   under.cex, under.font, under.ygap, ygap,
                   split.shadow.col, split.shadow.offset, bg,
                   small.underspace, split.strwidth, split.strheight)

        if(!draw.split.shadows1)
            draw.labs1(split.labs, split.boxes, split.yspace,
                split.cex * cex, split.font, split.family, split.col,
                draw.split.shadows1, draw.split.shadows, split.shadow.col, split.round)

        if(is.fancy(type)) { # right hand boxes
            right.split.boxes <-
                draw.boxes("right", draw.split.shadows1,
                   right.split.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   split.cex * cex, split.font, split.family, split.adj, split.yshift,
                   split.box.col, split.border.col,
                   split.lty, split.lwd,
                   split.space, split.yspace,
                   split.round * split.strheight,
                   under.cex, under.font, under.ygap, ygap,
                   split.shadow.col, split.shadow.offset, bg)

            if(!draw.split.shadows1)
                draw.labs1(right.split.labs, right.split.boxes, split.yspace,
                    split.cex * cex, split.font, split.family, split.col,
                    draw.split.shadows1, draw.split.shadows, split.shadow.col, split.round)
        }
        node.boxes <- draw.boxes("default", draw.shadows1,
                   node.labs, node.xy, xlim, ylim,
                   nodes, branch,
                   Margin, xflip, yflip, main, sub,
                   col.main, cex.main, col.sub, cex.sub,
                   cex, font, family, adj, yshift,
                   box.col, border.col,
                   lty, lwd,
                   space, yspace,
                   round * strheight,
                   under.cex, under.font, under.ygap, ygap,
                   shadow.col, shadow.offset, bg)

        draw.labs1(node.labs, node.boxes, yspace,
            cex, font, family, col,
            draw.shadows1, draw.shadows, shadow.col, round)

        if(yesno && !is.fancy(type) && !snip) # draw "yes" and "no" at root?
            rpart.plot:::draw.yes.no(yesno, yes.text, no.text,
                    type, draw.shadows1,
                    xflip, left, branch, xlim, ylim, node.xy, lwd,
                    yesno.yshift,
                    split.boxes, split.cex * cex, split.box.col, split.border.col,
                    split.shadow.col, split.shadow.offset,
                    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
                    nn.border.col, nn.lty, nn.round, bg)

        if(nn || ni)
            rpart.plot:::draw.node.numbers(nn, ni, draw.shadows1, type, branch,
                    Margin, xflip, yflip, cex,
                    main, sub, col.main, cex.main, col.sub, cex.sub,
                    xlim, ylim, node.xy, is.leaf, nodes,
                    node.labs, font,  family, box.col, border.col, shadow.col,
                    under.cex, under.font, under.ygap, ygap,
                    split.labs, split.cex * cex, split.font, split.family, split.box.col,
                    split.border.col, split.shadow.col,
                    nn.cex, nn.font, nn.family, nn.col, nn.box.col,
                    nn.border.col, nn.lty, nn.lwd, nn.round,
                    split.adj, split.space, split.yspace, split.yshift,
                    yshift, adj, space, yspace, shadow.offset,
                    nn.adj, nn.yshift, nn.space, nn.yspace, bg)

        list(node.boxes=node.boxes, split.boxes=split.boxes)
    }
    #--- prp starts here ---

    if(!inherits(x, "rpart"))
        stop("Not an rpart object")
    obj <- x

    # Process dots args.  The call to eval.parent is necessary to evaluate the
    # call to say "c" when user does something like "xlim=c(0,2)". Note
    # also that we allow abbreviations by using say "dots$fo" instead of "dots$font".
    # TODO Is there a better way? This approach is fragile, we have to be
    #       extremely careful that abbreviation doesn't alias with other args.
    # TODO If you specify col.main=2 it affects both main and the color of the nodes.

    old.warnPartialMatchDollar <- getOption("warnPartialMatchDollar")
    if(rpart.plot:::is.boolean(old.warnPartialMatchDollar)) # prevents problem when old value is NULL
        on.exit(options(warnPartialMatchDollar=old.warnPartialMatchDollar))
    options(warnPartialMatchDollar=FALSE)

    dots <- match.call(expand.dots=FALSE)$...
    check.dots(dots)
    adj      <- eval.parent(dots$adj);   if(is.null(adj))      adj      <- par("adj")
    cex.main <- eval.parent(dots$cex.m)
    cex.sub  <- eval.parent(dots$cex.s)
    col      <- eval.parent(dots$col);   if(is.null(col))      col      <- par("col")
    col.main <- eval.parent(dots$col.m); if(is.null(col.main)) col.main <- par("col.main")
    col.sub  <- eval.parent(dots$col.s); if(is.null(col.sub))  col.sub  <- par("col.sub")
    family   <- eval.parent(dots$fam);   if(is.null(family))   family   <- par("family")
    font     <- eval.parent(dots$fo);    if(is.null(font))     font     <- par("font")
    lty      <- eval.parent(dots$lt);    if(is.null(lty))      lty      <- par("lty")
    lwd      <- eval.parent(dots$lw);    if(is.null(lwd))      lwd      <- par("lwd")
    main     <- eval.parent(dots$mai)
    mar      <- eval.parent(dots$mar)
    sub      <- eval.parent(dots$sub);
    xlim     <- eval.parent(dots$xl)
    xpd      <- eval.parent(dots$xp)
    ylim     <- eval.parent(dots$yl)

    if(rpart.plot:::is.boolean(old.warnPartialMatchDollar))
        options(warnPartialMatchDollar=old.warnPartialMatchDollar)

    if(is.null(under.col))   under.col   <- col
    if(is.null(border.col))  border.col  <- col
    if(is.null(branch.lwd))  branch.lwd  <- lwd
    if(is.null(split.lwd))   split.lwd   <- lwd
    if(is.null(nn.lwd))      nn.lwd      <- lwd
    if(is.null(split.adj))   split.adj   <- adj

    class.stats <- NULL
    #if(obj$method == "class" || is.class.response(obj))
    #    class.stats <- get.class.stats(obj)

    # The idea with the following  argument checking is to catch user
    # errors here where possible before they cause an obscure message
    # later on.  But it is impossible to be exhaustive.

    trace <- as.numeric(rpart.plot:::check.numeric.scalar(trace, logical.ok=TRUE))
    type <- rpart.plot:::check.integer.scalar(type, logical.ok=FALSE)
    if(type < TYPE0.default || type > TYPE5.varname.in.node)
        stop("type must be ", TYPE0.default, "...",
              TYPE5.varname.in.node, ", you have type=", type)

    under <- rpart.plot:::check.boolean(under)
    clip.left.labs[1] <- rpart.plot:::check.boolean(clip.left.labs[1])
    clip.right.labs[1] <- rpart.plot:::check.boolean(clip.right.labs[1])
    nn <- rpart.plot:::check.boolean(nn)
    ni <- rpart.plot:::check.boolean(ni)
    stopifnot((is.numeric(yesno) || is.logical(yesno)) &&
              length(yesno) == 1 && floor(yesno) == yesno)
    if(yesno < 0 || yesno > 2)
        stop("yesno must be 0, 1, or 2.  You have yesno=", yesno)
    stopifnot(is.character(yes.text) && length(yes.text) == 1)
    stopifnot(is.character(no.text) && length(no.text) == 1)
    fallen.leaves <-rpart.plot:::check.boolean(fallen.leaves)
    clip.facs <-rpart.plot:::check.boolean(clip.facs)
    uniform <-rpart.plot:::check.boolean(uniform)
    left <-rpart.plot:::check.boolean(left)
    xflip <-rpart.plot:::check.boolean(xflip)
    yflip <-rpart.plot:::check.boolean(yflip)
    do.par <-rpart.plot:::check.boolean(do.par)
    snip <-rpart.plot:::check.boolean(snip)
    if(snip) {
        branch.col = "black"
        branch.lty = 1
    }
    compress <-rpart.plot:::check.boolean(compress)
    ycompress <-rpart.plot:::check.boolean(ycompress)
    xcompact <-rpart.plot:::check.boolean(xcompact)
    ycompact <-rpart.plot:::check.boolean(ycompact)
    add.labs <-rpart.plot:::check.boolean(add.labs)
    boxes.include.gap <-rpart.plot:::check.boolean(boxes.include.gap)
    stopifnot(all(split.round >= 0))
    stopifnot(all(nn.round >= 0))
    stopifnot(tweak > 0 && tweak <= 10) # upper limit is arb
    stopifnot(max.auto.cex >= 1)
    stopifnot(min(shift.amounts) >= 0 && max(shift.amounts) <= 10) # 10 is arb
    stopifnot(xcompact.ratio > 0 && xcompact.ratio <= 2) # 2 is arb, max useful is probably 1
    stopifnot(min.auto.cex >= 0 && min.auto.cex <= 1)
    stopifnot(branch >= 0 && branch <= 1)
    if(!is.null(snip.fun))
        rpart.plot:::check.func.args(snip.fun, "snip.fun", function(tree) NULL)
    if(length(family) != 1 || length(split.family) != 1 || length(nn.family) != 1)
        stop("prp: family argument must be length 1 (family cannot be vectorized)")
    digits   <- rpart.plot:::process.digits.arg(digits)
    varlen   <- rpart.plot:::check.integer.scalar(varlen, logical.ok=FALSE)
    faclen   <- rpart.plot:::check.integer.scalar(faclen, logical.ok=FALSE)
    roundint <- rpart.plot:::check.boolean(roundint)

    bg <- get.bg() # never returns NA or 0
    border.col       <- set.zero.to.bg(border.col,       bg)
    shadow.col       <- set.zero.to.bg(shadow.col,       bg)
    under.col        <- set.zero.to.bg(under.col,        bg)
    split.col        <- set.zero.to.bg(split.col,        bg)
    split.box.col    <- set.zero.to.bg(split.box.col,    bg)
    split.shadow.col <- set.zero.to.bg(split.shadow.col, bg)
    nn.col           <- set.zero.to.bg(nn.col,           bg)
    nn.box.col       <- set.zero.to.bg(nn.box.col,       bg)
    nn.border.col    <- set.zero.to.bg(nn.border.col,    bg)

    # We use parent.frame() here -- we can't use attr(,".Environment") from
    # terms(obj) because that gives the wrong environment if rpart called
    # from within a function.
    if(is.null(obj$varinfo)) # this "if" prevents duplicate warnings
        obj$varinfo <- rpart.plot:::get.modelframe.info(obj, roundint, trace,
                                           parent.frame(), "prp")
    if(!rpart.plot:::is.na.or.zero(branch.type)) {
        branch <- if(branch > .5) 1 else 0
        ycompact <- FALSE # want branches to be as vertical as possible
    }
    auto.cex <- FALSE
    if(is.null(cex)) {
        auto.cex <- TRUE    # automatically calculate cex
        cex <- 1
    }
    if(is.null(split.cex))
        split.cex <- 1
    if(fallen.leaves)
        compress <- FALSE
    if(!is.null(obj$frame$splits))
        stop("Old-style rpart object?  (frame$splits is NULL)")
    frame <- obj$frame
    is.leaf <- rpart.plot:::is.leaf(frame)
    nodes <- as.numeric(row.names(frame))

    if(is.auto(extra, n=1))
        extra <- get.default.extra(obj, class.stats)

    node.fun.name <- deparse(substitute(node.fun))
    node.labs <- internal.node.labs(obj, node.fun, node.fun.name, type, extra,
                                    under, xsep, digits, varlen,
                                    prefix, suffix, class.stats, under.percent)

    # handle the box.col and box.palette arguments possibly specified by the user
    ret <- rpart.plot:::handle.box.palette.args(obj, trace, box.col, box.palette,
                                   pal.thresh, pal.node.fun,
                                   node.fun.name, class.stats, node.labs)
    box.col     <- ret$box.col     # box.palette (if specified) converted to box.col
    box.palette <- ret$box.palette # expanded box.palette
    box.col <- rpart.plot:::recycle(box.col, node.labs)
    if(type == TYPE5.varname.in.node) {
        box.col[!is.leaf] <- split.box.col
        col <- rpart.plot:::recycle(col, is.leaf)
        col[!is.leaf] <- split.col
        border.col <- rpart.plot:::recycle(border.col, is.leaf)
        border.col[!is.leaf] <-
            if(rpart.plot:::is.specified(split.border.col)) split.border.col else 1
    }
    split.labs <- split.labs.wrapper(obj, split.fun,
                deparse(substitute(split.fun)),
                split.prefix, split.suffix,
                right.split.prefix, right.split.suffix,
                type, clip.facs, clip.left.labs, clip.right.labs, xflip,
                digits, varlen, faclen, roundint, trace,
                facsep, eq, lt, ge)

    if(do.par) {
        # Make the side margins small.
        # Retain the top edge for the main title but only if necessary.
        # Likewise the bottom edge for the subtitle.
        # Note that family may change in my.strheight and init.plot, so we on.exit it here.
        init.plot(1, 1, Margin, xflip, yflip, main, sub,
                  col.main, cex.main, col.sub, cex.sub)
        par <- par("mar", "xpd", "family")
        on.exit(par(par))
        if(is.null(mar)) { # user did not explictly set mar when invoking prp?
            mar <- par$mar
            if(is.null(sub))  mar[1] <- 1
            if(is.null(main)) mar[3] <- 1
            mar[2] <- mar[4] <- 1
        }
        if(is.null(xpd)) # user did not explicitly set xpd when invoking prp?
            xpd <- NA
        par(mar=mar, xpd=xpd)
        par(new=TRUE) # done par for now, start next plot on the same page
    }
    if(is.fancy(type)) {
        right.split.labs <- split.labs[match(2 * nodes+1, nodes)]
        split.labs <- split.labs[match(2 * nodes, nodes)]
        if(!left) # TODO msg uses hard coded TYPE3.fancy, TYPE4.fancy.all, TYPE5.varname.in.node
            stop("left=FALSE is not yet supported with type=3 or 4 or 5")
    } else {
        if(left != xflip)   # default, set right labs to NA
            split.labs <- split.labs[match(2 * nodes, nodes)]
        else                # set left labs to NA
            split.labs <- split.labs[match(2 * nodes+1, nodes)]
    }
    draw.shadows <- !rpart.plot:::is.invisible(shadow.col, bg)
    draw.split.shadows <- !rpart.plot:::is.invisible(split.shadow.col, bg)

    # Recycle stuff that doesn't get recyled automatically.  It's more efficient
    # torpart.plot:::recycle it here once rather than over and over in get.boxes etc.
    adj                 <-rpart.plot:::recycle(adj, nodes)
    space               <-rpart.plot:::recycle(space, nodes)
    yspace              <-rpart.plot:::recycle(yspace, nodes)
    shadow.offset       <-rpart.plot:::recycle(shadow.offset, nodes)
    under.cex           <-rpart.plot:::recycle(under.cex, nodes)
    under.ygap          <-rpart.plot:::recycle(under.ygap, nodes)
    split.cex           <-rpart.plot:::recycle(split.cex, nodes)
    split.adj           <-rpart.plot:::recycle(adj, nodes)
    split.space         <-rpart.plot:::recycle(split.space, nodes)
    split.yspace        <-rpart.plot:::recycle(split.yspace, nodes)
    split.shadow.offset <-rpart.plot:::recycle(split.shadow.offset, nodes)
    nn.adj              <-rpart.plot:::recycle(nn.adj, nodes)
    nn.space            <-rpart.plot:::recycle(nn.space, nodes)
    nn.yspace           <-rpart.plot:::recycle(nn.yspace, nodes)

    ret <- get.yshift(type, nodes, is.leaf,
                      cex, node.labs, yshift, yspace, under.cex,
                      split.labs, split.cex, split.yshift, split.yspace, ygap)
    yshift       <- ret$yshift
    split.yshift <- ret$split.yshift

    if(yesno == 2 && !is.fancy(type))
        split.labs <- ifelse(split.labs == "NA",
                             "NA", paste(yes.text, split.labs, no.text))

    layout <- rpart.plot:::get.layout(obj, type, nn, yesno, fallen.leaves, branch,
        uniform, Margin, cex, auto.cex, compress, ycompress,
        trace, main, sub,
        node.labs, font, family, box.col, border.col,
        under.font, under.cex,
        split.labs, right.split.labs, split.cex, split.font, split.family,
        split.box.col, split.border.col,
        nspace, minbranch, adj, yshift, space, yspace,
        split.adj, split.yshift, split.space, split.yspace,
        gap, ygap, under.ygap, xcompact, ycompact, xcompact.ratio,  min.inter.height,
        max.auto.cex, min.auto.cex, ycompress.cex, accept.cex,
        shift.amounts, Fallen.yspace, bg)

    if(yesno == 2 && !is.fancy(type)) # remove prepended yes.text and no.text?
        split.labs <-
            ifelse(split.labs == "NA",
                     "NA",
                     substr(split.labs, nchar(yes.text)+1, nchar(split.labs)-nchar(no.text)-1))
    cex <- layout$cex
    gap <- layout$gap
    ygap <- layout$ygap
    # we use pmax here so there is always a little space even if it causes overlapping
    space <- pmax(.25, layout$node.space)
    yspace <- pmax(.25, layout$node.yspace)
    if(is.null(xlim))
        xlim <- layout$xlim
    stopifnot(is.numeric(xlim) && length(xlim) == 2)
    if(is.null(ylim))
        ylim <- layout$ylim
    stopifnot(is.numeric(ylim) && length(ylim) == 2)
    split.yshift <- layout$split.yshift
    if(trace > 0) {
        tweak.msg <- if(tweak == 1) "" else rpart.plot:::sprint(" (before applying tweak %g)", tweak)
        #printf("cex %.3g%s   xlim c(%.3g, %.3g)   ylim c(%.3g, %.3g)\n",
         #      cex[1], tweak.msg, xlim[1], xlim[2], ylim[1], ylim[2])
        print("cex %.3g%s   xlim c(%.3g, %.3g)   ylim c(%.3g, %.3g)\n",
               cex[1], tweak.msg, xlim[1], xlim[2], ylim[1], ylim[2])
    }
    if(!auto.cex && tweak != 1)
        warning("cex and tweak both specified, applying both")
    cex <- tweak * cex
    all.cex <- merge1(cex, split.cex * cex)

    split.lwd  <-rpart.plot:::recycle(cex * split.lwd,  nodes)
    branch.lwd <-rpart.plot:::recycle(cex * branch.lwd, nodes)
    nn.lwd     <-rpart.plot:::recycle(cex * nn.lwd,     nodes)
    # do this last because split.lwd etc. above use lwd as the default
    lwd <-rpart.plot:::recycle(cex * lwd, nodes)

    node.xy <- layout$node.xy
    init.plot(xlim, ylim, Margin, xflip, yflip, main, sub,
              col.main, cex.main, col.sub, cex.sub,
              fam.main=fam.main, cex=cex[1], trace=trace, hide.title=FALSE)
    split.strwidth  <- rpart.plot:::my.strwidth("M", split.cex * cex, split.font, split.family)
    strheight <- rpart.plot:::my.strheight("M", cex, font, family)
    split.strheight <- rpart.plot:::my.strheight("M", split.cex * cex, split.font, split.family)
    node.boxes <- split.boxes <- NA
    if(add.labs) {
        if(is.null(round))
            round <- max(1, 2 * min(space, yspace))
        stopifnot(all(round >= 0))
        round <-rpart.plot:::recycle(round, nodes)
        if(is.null(leaf.round))
            leaf.round <- round
        stopifnot(all(leaf.round >= 0))
        leaf.round <-rpart.plot:::recycle(leaf.round, nodes)
        round[is.leaf] <- leaf.round[is.leaf]
        if(type == TYPE5.varname.in.node) # want leaf boxes to look different to interior boxes
            round[!is.leaf] <- 0           # TODO would be nice if this was more flexible
        # draw shadows first, if any, so boxes and lines are over shadows
        if(draw.shadows || draw.split.shadows)
            draw.labs(draw.shadows, draw.split.shadows)
    }
    # draw branch lines now so they are over the shadows, under the boxes
    branch.xy <- rpart.plot:::draw.branches(obj, branch.type, branch.col,
                    branch.lty, branch.lwd,  branch.fill, branch.tweak,
                    node.labs, split.labs, node.xy, strheight,
                    type, branch, xflip, yflip, Margin, space, yspace,
                    cex, font, family, adj, box.col, border.col,
                    split.cex, split.font, split.family, split.adj, split.yshift,
                    split.box.col, split.border.col, split.space, split.yspace,
                    main, sub, col.main, cex.main, col.sub, cex.sub,
                    xlim, ylim, yshift, ygap, bg,
                    min.branch.width)
    if(add.labs) {
        ret <- draw.labs(FALSE, FALSE)
        node.boxes  <- ret$node.boxes
        split.boxes <- ret$split.boxes
    }
    snipped.nodes <- NULL
    #if(snip) {
    #    ret <- do.snip(obj, nodes, split.labs, node.xy, branch.xy,
    #                   branch.lwd, xlim, ylim, digits, snip.fun, cex)
    #    obj <- ret$obj
    #    snipped.nodes <- ret$snipped.nodes
    #}
    ret <- list(obj=obj, snipped.nodes=snipped.nodes,
                xlim=xlim, ylim=ylim,
                x=node.xy$x, y=node.xy$y,
                branch.x=branch.xy$x, branch.y=branch.xy$y,
                labs=node.labs, cex=cex, boxes=node.boxes,
                split.labs="", split.cex=split.cex, split.box=split.boxes)

    rpart.plot:::possible.palette.legend(ret, class.stats, box.col, box.palette,
                            legend.x, legend.y, legend.cex, tweak, trace)

    invisible(ret)
}

#' Another rpart.plot function
#' @keywords internal
#' @noRd
#' @importFrom graphics axis grid par plot rect text
init.plot <- function(x, y,
                      Margin, xflip, yflip, main, sub,
                      col.main, cex.main, col.sub, cex.sub,
                      fam.main="", cex=1, trace=0, hide.title=TRUE)
{
    if(length(x) == 1)
        x <- c(0, x)
    if(length(y) == 1)
        y <- c(0, y)
    xlim <- range(x) + diff(range(x)) * c(-Margin, Margin)
    if(xflip)
        xlim <- rev(xlim)
    ylim <- range(y) + diff(range(y)) * c(-Margin, Margin)
    if(yflip)
        ylim <- rev(ylim)
    old.family <- NA
    if(hide.title) {
        # need to plot main and sub to get screen layout that accounts
        # for their allotted space, but want them to be invisible
        col.main <- col.sub <- 0
    } else {
        if(is.null(cex.main))
            cex.main <- max(.8, min(1.5 * cex, 1.2))
        if(is.null(cex.sub))
            cex.sub <- max(.8, min(1.5 * cex, 1.2))
        if(!is.null(sub)) # hack to make subtitle visible with our do.par()
            sub <- paste(sub, "\n")
        if(!identical(fam.main, ""))
            old.family <- par(family=fam.main) # on.exit for this already set up in prp()
    }
    plot(0, 0, xlim=xlim, ylim=ylim, type="n", axes=FALSE, xlab="", ylab="",
         main=main, sub=sub,
         col.main=col.main, cex.main=cex.main, col.sub=col.sub, cex.sub=cex.sub)
    if(!is.na(old.family))
        par(family=old.family)
    if(trace >= 2) { # draw the grid and region boxes
        col <- "palegreen"
        # set xpd so grid lines stay in our region
        old.xpd <- par(xpd=FALSE)
        on.exit(par(xpd=old.xpd))
        grid(col=col, lwd=cex)
        axis(1, col=col, col.axis=col)
        axis(2, col=col, col.axis=col)
        rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border=col, lty=1, lwd=cex)
        text((xlim[1] + xlim[2]) / 2, ylim[2], "xlim ylim", col=col)
        usr <- par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], col=NA, border=col, lty=1, lwd=cex)
        text((usr[1] + usr[2]) / 2, usr[4], "usr", col=col, xpd=NA)
    }
}

#' This is the type of thing I think I can probably delete.
#process.digits.arg <- function(digits)
#{
#    digits <- rpart.plot:::check.integer.scalar(digits, min=-22, max=22, logical.ok=FALSE)
#    if(digits == 0)
#        digits <- getOption("digits")
#    if(digits > 8)  # silently reduce digits because of verysmall in tweak.splits
#        digits <- 8
#    else if(digits < -8)
#        digits <- -8
#    digits
#}


#' @keywords internal
#' @noRd
get.default.extra <- function(obj, class.stats)
{
    #if(obj$method == "class" || is.class.response(obj)) {
     #   if(class.stats$nlev > 2)
       #     104 # multiclass response
      #  else
       #     106 # binomial model (two class response)
    #} else
  #if(obj$method == "poisson" || obj$method == "exp")
  #      101
  #  else
        100
}

#' @keywords internal
#' @noRd
get.yshift <- function(type, nodes, is.leaf,
                       cex, node.labs, yshift, yspace, under.cex,
                       split.labs, split.cex, split.yshift, split.yspace, ygap)
{
    # Return number of lines (separated by \n) in each lab in labs
    # A \n\n counts as one \n (it should really equal split.yshift)
    # TODO there must be a simpler way to do this
    get.nlines <- function(labs)
    {
        labs <- gsub("\n\n", "\n", labs) # replace \n\n with \n
        sapply(strsplit(labs, "\n"), function(x) length(x))
    }
    #--- get.yshift starts here ---
    yshift       <-rpart.plot:::recycle(yshift, nodes)
    split.yshift <-rpart.plot:::recycle(split.yshift, nodes)

    # always want _some_ ygap else cex has to be very small to prevent overplotting
    if(length(ygap) == 0)
        ygap <- .5

    # We use the number of lines to estimate the vert space taken by the label
    # We discount the first line of the inbox text.

    node.sep.labs <- separate.labs(node.labs)
    node.nlines  <- get.nlines(node.sep.labs$in.box) - 1 +
                    under.cex * get.nlines(node.sep.labs$under.box)

    split.sep.labs <- separate.labs(split.labs)
    split.nlines  <- get.nlines(split.sep.labs$in.box) - 1 +
                    under.cex * get.nlines(split.sep.labs$under.box)

    # Note that get.boxes uses split.yshift with split.cex, but need yshift
    # w.r.t. node.cex --- hence the node.to.split conversion ratio.

    ratio <- cex / split.cex

    # The following calculations must match the y1 and y2 calculations in get.boxes
    # We have to calculate the position of the _center_ of the boxes.
    # TODO investigate .9 etc. below (empirically determined) --- space for newline?

    node.shift  <- -yshift + .6 + .9 * node.nlines + .5 * yspace
    split.shift <- .5 * yspace + .6 + .9 * split.nlines

    if(is.fancy(type)) {
        # Want the node box on the node, and the top of left split box below
        # the bottom of node box, with a little vertical space.
        # The right split box position is calculated in get.boxes, once we
        # have a known cex and calculate the absolute box sizes.
        # The node box, left split and right split get treated as three
        # separate boxes in get.layout.

        split.yshift <- split.yshift -
                        ratio * node.shift - max(ygap, .25) - split.shift

    } else if(type == TYPE1.all) {
        # Want the split box on the node, and the top of the node
        # box just above the bottom of the split box (a slight overlap).
        # These get combined and treated as one large box in get.layout.

        new.yshift <- -node.shift - split.shift / ratio

        # following needed when user explicitly sets yshift
        split.yshift <- split.yshift + yshift
        # we only want to shift node labels with a split label sitting on them
        yshift[!is.leaf] <- new.yshift[!is.leaf]
        # except that we do move leaf nodes down a bit to roughly match their split brothers
        yshift[is.leaf] <- yshift[is.leaf] - node.nlines[is.leaf]

    } else if(type == TYPE2.all.under) {
        # Want the node box on the node, and the top of the split
        # box just below the bottom of the node box.
        # These get combined and treated as one large box in get.layout.
        split.yshift <- split.yshift -
                        ratio * node.shift - split.shift
    }
    list(yshift=yshift, split.yshift=split.yshift)
}

#' @keywords internal
#' @noRd
get.node.coords <- function(obj, uniform, branch, compress,
                            nspace, minbranch, fallen.leaves, Fallen.yspace)
{
    if(NROW(obj$frame) <= 1)    # tree is just a root?
        compress <- FALSE       # prevent rpartco from issuing an error
    # initialize nspace for rpartco()
    if(!compress)
        nspace <- -1    # magic value for rpartco meaning no compression
    if(is.null(nspace))
        nspace <- branch
    xy <- rpart.plot:::my.rpartco(obj, uniform, nspace, minbranch)
    x <- xy$x
    y <- xy$y

    # scale x to 0 ... 1
    x <- x - min(x)
    if(length(x) > 1)   # needed when tree is just a root
        x <- x  / max(x)

    # scale y to 0 ... 1 so top node is always at 1, bottom leaves at 0
    y <- y - min(y)
    if(length(y) > 1)
        y <- y  / max(y)

    if(fallen.leaves) {
        is.leaf <- rpart.plot:::is.leaf(obj$frame)
        y[is.leaf] <- 0
        # Make more space above leaves by shifting all other nodes up to make to
        # more likely that we will be able to shift fallen leaves for more space.
        # It usually looks a little better too.
        y[!is.leaf] <- (y[!is.leaf] + Fallen.yspace) / (1 + Fallen.yspace)
    }
    list(x=x, y=y)
}
# Get the box coords, a row for each box.
# Use do.init.plot=FALSE if want to use char sizes etc. of the existing plot.


#' @keywords internal
#' @noRd
get.boxes <- function(boxtype,  # one of "default", "left", "right", "undersplit"
    labs, x, y, xlim, ylim, nodes, branch,
    Margin, xflip, yflip, main, sub, col.main, cex.main, col.sub,  cex.sub,
    cex, font, family, adj, yshift, box.col, border.col, space, yspace,
    ygap, bg, do.init.plot=TRUE,
    box.around.all.text=TRUE)   # else box only around "in box" text i.e. text before \n\n
{                               # TRUE when figuring out box spacing, FALSE when drawing boxes
    if(do.init.plot)
        init.plot(xlim, ylim, Margin, xflip, yflip, main, sub,
                  col.main, cex.main, col.sub, cex.sub)

    # to minimize blanking out parts of the branch lines, we want only a
    # small white space around letters when fancy and non-visible boxes
    if((boxtype == "left" || boxtype == "right") &&
            rpart.plot:::is.box.invisible(box.col, border.col, bg))
        space <-rpart.plot:::recycle(min(.2, space), labs)

    # sanity check that variables are already expanded correctly for recycling
    stopifnot(length(adj)    == length(labs) &&
              length(yshift) == length(labs) &&
              length(space)  == length(labs) &&
              length(yspace) == length(labs))

    # TODO simplistic approach for now, assumes \n approx equal to under.ygap
    stripped.labs <- gsub("\n\n", "\n", labs) # replace \n\n with \n
    sep.labs <- separate.labs(labs)
    in.box.labs <- sep.labs$in.box

    height1         <-rpart.plot:::recycle(rpart.plot:::my.strheight("M", cex, font, family), labs)
    stripped.height <- rpart.plot:::my.strheight(stripped.labs, cex, font, family)
    width1          <-rpart.plot:::recycle(rpart.plot:::my.strwidth("M", cex, font, family), labs)
    width <- if(box.around.all.text)
                rpart.plot:::my.strwidth(labs,        cex, font, family)
             else
                rpart.plot:::my.strwidth(in.box.labs, cex, font, family)
    if(do.init.plot)
        par(new=TRUE) # done strwidth and strheight, ensure we stay on same page
    xy.to.calc.xshift <- list(x=x, y=y) # grab these before they change

    # The text function (which draws the labels elsewhere), centers the labels
    # vertically on the point, a fact that is used to calculate box positions below.
    # The following calculations must match the calculations in get.yshift.

    x[is.na(labs)] <- y[is.na(labs)] <- NA
    x1 <- x - adj     * width - space/2 * width1                      # left edge
    x2 <- x + (1-adj) * width + space/2 * width1                      # right edge
    y2 <- y + yshift * height1 + stripped.height/2 + yspace * height1 # top of box
    if(box.around.all.text)
        y1 <- y2 - stripped.height - yspace * height1
    else {
        in.box.height <- rpart.plot:::my.strheight(in.box.labs, cex, font, family)
        y1 <- y2 - in.box.height - yspace * height1
    }
    xshift <- 0 # stays at zero unless left or right split
    if(length(x1) > 1) { # tree is not just a root?
        if(boxtype == "undersplit") {
            # splits are under the node boxes,
            # force a little extra space under the split box so branch line is visible
            stopifnot(box.around.all.text) # should use this only in get.layout
            y2 <- y2 - height1
        } else if(boxtype == "left") {
            child <- match(2 * nodes, nodes)        # left child
        } else if(boxtype == "right") {
            child <- match(2 * nodes + 1, nodes)    # right child
            # lower the right splits relative to the left splits
            box.heights <- y2 - y1
            # May 2018: Raised right hand labels slightly when actually drawing
            # because they sometimes were too close to the node box beneath them.
            # We only do it when actually drawing otherwise the layout engine
            # code elsewhere would have to be adjusted.
            adjust <-
                if(box.around.all.text) # figuring out spacing?
                    box.heights + (max(ygap, .4) * height1) # (original code)
                else                    # actually drawing
                    box.heights
            y1 <- y1 - adjust
            y2 <- y2 - adjust
            yshift <- yshift - min(box.heights / height1, na.rm=TRUE)
        }
        if(boxtype == "left" || boxtype == "right") {
            # adjust x coords so labels are centered on the branch lines
            branch.xy <- rpart.plot:::get.branches(xy.to.calc.xshift$x, xy.to.calc.xshift$y, nodes, branch)
            x <- branch.xy$x[, child]
            y <- branch.xy$y[, child]
            xshift <- (x[2, ] - x[3, ]) + # 1.3 below to exaggerate the separation, looks better
                      1.3 * yshift * height1 * (x[2, ] - x[1, ]) / (y[2, ] - y[1, ])
        }
    }
    list(x1=x1 + xshift, y1=y1, x2=x2 + xshift, y2=y2)
}

#' @keywords internal
#' @noRd
draw.boxes <- function(fancy.style, draw.shadow, labs, xy,
                       xlim, ylim, nodes, branch,
                       Margin, xflip, yflip, main, sub,
                       col.main, cex.main, col.sub, cex.sub,
                       cex, font, family, adj, yshift,
                       box.col, border.col,
                       lty, lwd, space, yspace, r,
                       under.cex, under.font, under.ygap, ygap,
                       shadow.col, shadow.offset, bg,
                       small.underspace=FALSE, split.strwidth=0, split.strheight=0)
{
    box <- get.boxes(fancy.style, labs, xy$x, xy$y, xlim, ylim, nodes, branch,
                     Margin, xflip, yflip, main, sub,
                     col.main, cex.main, col.sub, cex.sub,
                     cex, font, family, adj,
                     yshift, box.col, border.col, space, yspace,
                     ygap, bg,
                     do.init.plot=FALSE,
                     box.around.all.text=FALSE)

    new.box <- box
    if(small.underspace) {
        # Splits are under the node boxes.  Reduce sides and bottom of box slightly so
        # just a little white space below the split box so branch line is more visible.
        add.space  <- pmin(space, .6)
        add.yspace <- pmin(yspace, .8)
        new.box$x1 <- new.box$x1 + (space - add.space)   * split.strwidth
        new.box$x2 <- new.box$x2 - (space - add.space)   * split.strwidth
        new.box$y1 <- new.box$y1 + (yspace - add.yspace) * split.strheight
    }
    if(!draw.shadow)
        rpart.plot:::rounded.rect(new.box$x1, new.box$y1, new.box$x2, new.box$y2,
                     xlim, ylim, r, box.col, border.col, lty, lwd)
    else if(!rpart.plot:::is.invisible(shadow.col, bg))
        rpart.plot:::draw.shadow(new.box$x1, new.box$y1, new.box$x2, new.box$y2,
                    xlim, ylim, r, shadow.col, shadow.offset)
    box
}

#' @keywords internal
#' @noRd
# Set bg to the background color or "white" if transparent.
# The idea is that we want a color that is opaque but matches background.
get.bg <- function()
{
    bg <- par("bg")
    if(bg[1] == "transparent" || # TODO par("bg") incorrectly(?) returns transparent with mfrow
       bg[1] == 0 || is.na(bg[1])) { # probably unnecessary
        bg <- "white"
    }
    bg
}

#' @keywords internal
#' @noRd
set.zero.to.bg <- function(col, bg) # set elems of col that are 0 or NA to bg
{
    if(is.null(col))
        col <- bg
    else
        col[which(col == 0) | is.na(col)] <- bg
    col
}


#' @keywords internal
#' @noRd
# true if x == "auto" or "-auto", ignoring case, partial match to n characters
is.auto <- function(x, n=2)
{
    is.character(x) &&
    length(x) >= 1  &&
    if(n == 1) # only one character is needed to disambiguate from e.g. extra=1
        grepl("^a",     substr(x[1], 1, 1), ignore.case=TRUE)
    else       # two characters needed to disambiguate from e.g. "aliceblue"
        (grepl("^au",   substr(x[1], 1, 2), ignore.case=TRUE) ||
         grepl("^\\-a", substr(x[1], 1, 2), ignore.case=TRUE))
}

#' @keywords internal
#' @noRd
is.fancy <- function(type)
{
    type == TYPE3.fancy.no.all ||
    type == TYPE4.fancy.all    ||
    type == TYPE5.varname.in.node
}

# text before \n\n goes in the box
# text after \n\n if any goes under the box
#' @keywords internal
#' @noRd
separate.labs <- function(labs) {
    labs <- strsplit(labs, "\n\n")
    list(in.box    = sapply(labs, function(x) x[1]), under.box = sapply(labs, function(x) paste(x[-1], collapse="\n")))
}

#' @keywords internal
#' @noRd
get.box.centers <- function(box) {list(x=(box$x1 + box$x2)/2, y=(box$y1 + box$y2)/2)}

#' @keywords internal
#' @noRd
# split.labs.R: functions for generating split.labels
split.labs.wrapper <- function(x, split.fun, split.fun.name,
                               split.prefix, split.suffix,
                               right.split.prefix, right.split.suffix,
                               type, clip.facs,
                               clip.left.labs, clip.right.labs, xflip,
                               digits, varlen, faclen, roundint, trace,
                               facsep, eq, lt, ge)
{
    logical.eq <- eq
    if(clip.facs)
        eq <- "|" # special value used as a flag

    labs <- internal.split.labs(x, type,
                                digits, varlen, faclen, roundint,
                                clip.facs, clip.left.labs, clip.right.labs, xflip,
                                trace,
                                facsep, eq, logical.eq, lt, ge,
                                split.prefix, right.split.prefix,
                                split.suffix, right.split.suffix)

    if(!is.null(split.fun)) { # call user's split.fun?
        rpart.plot:::check.func.args(split.fun, "split.fun",
                        function(x, labs, digits, varlen, faclen) NA)
        labs <- split.fun(x, labs, abs(digits), varlen, faclen)
    }
    # check returned labs because split.fun may be user supplied
    rpart.plot:::check.returned.labs(x, labs, split.fun.name)
    labs
}

#' @keywords internal
#' @noRd
# Modified version of labels.rpart.
# This uses format0 instead of formatg and has various other extensions.
internal.split.labs <- function(x, type,
                                digits, varlen, faclen, roundint,
                                clip.facs, clip.left.labs, clip.right.labs, xflip,
                                trace,
                                facsep, eq, logical.eq, lt, ge,
                                split.prefix, right.split.prefix,
                                split.suffix, right.split.suffix)
{
    frame <- x$frame
    if(nrow(frame) == 1) # special case, no splits?
        return("root")   # NOTE: return
    is.leaf <- rpart.plot:::is.leaf(frame)
    split.var.names <- frame$var[!is.leaf]  # variable names for the (primary) splits

    split.var.names <- as.character(split.var.names) # factor levels to character
    clip.left.labs  <-rpart.plot:::recycle(clip.left.labs,  split.var.names)
    clip.right.labs <-rpart.plot:::recycle(clip.right.labs, split.var.names)

    # isplit is the row index of the primary split in x$splits
    index <- cumsum(c(1, frame$ncompete + frame$nsurrogate + !is.leaf))
    isplit  <- index[c(!is.leaf, FALSE)]

    split <- get.lsplit.rsplit(x, isplit, split.var.names,
                               type, clip.left.labs, clip.right.labs, xflip,
                               digits, faclen, roundint, trace,
                               facsep, eq, logical.eq, lt, ge)

    # We now have something like this:
    #    split.var.names:   sex   age     pclass     sibsp   pclass
    #    split$lsplit:      mal   >=9.5   =2nd,3rd   >=2.5   =3rd
    #    split$rsplit:      fml   <9.5    =1st       <2.5    =1st,2nd

    if(clip.facs) {
        is.eq.l <- substr(split$lsplit, 1, 1) == "|" # special value used a flag
        is.eq.r <- substr(split$rsplit, 1, 1) == "|"
        split.var.names[is.eq.l | is.eq.r] <- ""                     # drop var name
        split$lsplit[is.eq.l] <- substring(split$lsplit[is.eq.l], 2) # drop leading |
        split$rsplit[is.eq.r] <- substring(split$rsplit[is.eq.r], 2)
    }
    paste.split.labs(frame,
                     split.var.names, split$lsplit, split$rsplit,
                     type, clip.facs,
                     clip.left.labs, clip.right.labs, xflip, varlen,
                     split.prefix, right.split.prefix,
                     split.suffix, right.split.suffix)
}

#' @keywords internal
#' @noRd
get.lsplit.rsplit <- function(x, isplit, split.var.names,
                              type, clip.left.labs, clip.right.labs, xflip,
                              digits, faclen, roundint, trace,
                              facsep, eq, logical.eq, lt, ge)
{
    frame <- x$frame
    is.leaf <- rpart.plot:::is.leaf(frame)
    splits <- rpart.plot:::tweak.splits(x, roundint, digits, trace)
    ncat  <- splits[isplit, "ncat"]
    lsplit <- rsplit <- character(length=length(isplit))
    is.con <- ncat <= 1             # continuous vars (a logical vector)
    if(any(is.con)) {               # any continuous vars?
        cut <- splits[isplit[is.con], "index"]
        formatted.cut <- rpart.plot:::format0(cut, digits)
        is.less.than <- ncat < 0
        lsplit[is.con] <- paste0(ifelse(is.less.than, lt, ge)[is.con], formatted.cut)
        rsplit[is.con] <- paste0(ifelse(is.less.than, ge, lt)[is.con], formatted.cut)
        # print logical and 01 predictors as "Survived = 1" or "Survived = 0"
        islogical <- x$varinfo$islogical[isplit]
        is01      <- x$varinfo$is01[isplit]
        if(!anyNA(islogical) && !anyNA(is01) && (any(islogical) || any(is01))) {
            isbool <- islogical | (roundint & is01)
            eq0 <- paste0(logical.eq, "0")
            eq1 <- paste0(logical.eq, "1")
            lsplit[isbool] <- paste0(ifelse(is.less.than, eq0, eq1)[isbool])
            rsplit[isbool] <- paste0(ifelse(is.less.than, eq1, eq0)[isbool])
        }
    }
    is.cat <- ncat > 1              # categorical variables (a logical vector)
    if(any(is.cat)) {               # any categorical variables?
        # jrow is the row numbers of factors within lsplit and rsplit
        # cindex is the index on the "xlevels" list
        jrow <- seq_along(ncat)[is.cat]
        crow <- splits[isplit[is.cat], "index"] # row number in csplit
        xlevels <- attr(x, "xlevels")
        cindex <- match(split.var.names, names(xlevels))[is.cat]
        # decide if we must add a "=" prefix
        paste.left.eq  <- !is.fancy(type) | (if(xflip) !clip.right.labs else !clip.left.labs)
        paste.right.eq <- !is.fancy(type) | (if(xflip) !clip.left.labs  else !clip.right.labs)
        for(i in 1:length(jrow)) {
            node.xlevels <- rpart.plot:::my.abbreviate(xlevels[[cindex[i]]],
                                          faclen, one.is.special=TRUE)
            j <- jrow[i]
            splits <- x$csplit[crow[i],]
            # splits is 1=left 2=neither 3=right
            left  <- (1:length(splits))[splits == 1]
            right <- (1:length(splits))[splits == 3]
            collapse <- if(faclen==1) "" else facsep
            lsplit[j] <- paste(node.xlevels[left],   collapse=collapse)
            rsplit[j] <- paste0(node.xlevels[right], collapse=collapse)
            if(paste.left.eq[i])
                lsplit[j] <- paste0(eq, lsplit[j])
            if(paste.right.eq[i])
                rsplit[j] <- paste0(eq, rsplit[j])
        }
    }
    list(lsplit=lsplit, rsplit=rsplit)
}

# Paste the various constituents to create the split labels vector.
# On entry we have something like this:
#    split.var.names sex   age     pclass     sibsp   pclass
#    lsplit:         mal   >=9.5   =2nd,3rd   >=2.5   =3rd
#    rsplit:         fml   <9.5    =1st       <2.5    =1st,2nd
#' @keywords internal
#' @noRd
paste.split.labs <- function(frame, split.var.names, lsplit, rsplit,
                             type, clip.facs,
                             clip.left.labs, clip.right.labs, xflip, varlen,
                             split.prefix, right.split.prefix,
                             split.suffix, right.split.suffix)
{
    nodes <- as.numeric(row.names(frame))
    is.right <- nodes %% 2 == 1
    is.leaf <- rpart.plot:::is.leaf(frame)
    parent <- match(nodes %/% 2, nodes[!is.leaf])
    split.var.names <- rpart.plot:::my.abbreviate(split.var.names, varlen)
    left.names <- right.names <- split.var.names

    if(is.fancy(type)) {
        if(xflip)
            right.names[clip.left.labs] <- ""
        else
            left.names[clip.left.labs] <- ""
        if(xflip)
            left.names[clip.right.labs] <- ""
        else
            right.names[clip.right.labs] <- ""
    }
    if(is.null(right.split.prefix))
        right.split.prefix <- split.prefix
    if(is.null(right.split.suffix))
        right.split.suffix <- split.suffix

    # the heart of this function
    newline <- "\n\n"

    labs  <- paste0(split.prefix, left.names[parent], lsplit[parent], split.suffix,
                    newline, "pval", frame$pval)

    labs[is.right] <- paste0(right.split.prefix,
                             right.names[parent],
                             rsplit[parent],
                             right.split.suffix)[is.right]
    labs[1] <- "root" # was "NANA" because of above paste0
    labs
}




# node.labs.R: functions for generating labels
EX0                                 <- 0
EX1.NOBS                            <- 1
EX2.CLASS.RATE                      <- 2
EX3.MISCLASS.RATE                   <- 3
EX4.PROB.PER.CLASS                  <- 4
EX5.PROB.PER.CLASS.DONT             <- 5
EX6.PROB.2ND.CLASS                  <- 6
EX7.PROB.2ND.CLASS.DONT             <- 7
EX8.PROB.FITTED.CLASS               <- 8
EX9.PROB.ACROSS.ALL                 <- 9
EX10.PROB.ACROSS.ALL.2ND.CLASS      <- 10
EX11.PROB.ACROSS.ALL.2ND.CLASS.DONT <- 11


# call node.fun or obj$functions$text, and check its args and returned value.
# I actually modified this one as well!!!!
#' @keywords internal
#' @noRd
internal.node.labs <- function(x, node.fun, node.fun.name, type, extra,
                               under, xsep, digits, varlen,
                               prefix, suffix, class.stats, under.percent)
{
    stopifnot(is.numeric(extra) || is.logical(extra))
    stopifnot(length(extra) == 1)
    ex <- if(extra < 100) extra else extra - 100
    if(ex < 0 || floor(ex) != ex) {
        #extra.help()
        stop("extra=", extra, " is illegal")
    }
    stopifnot(extra >= 0)
    stopifnot(is.character(x$method) && length(x$method) == 1) # sanity check
    frame <- x$frame
    labs <-
        if(x$method == "anova")
            get.anova.labs(x, extra, under, digits, xsep, varlen, under.percent)
    else {
           stop("Unrecognized rpart object")
    }
    if(!is.null(node.fun)) {
        # call user's node.fun
        node.fun <- rpart.plot:::check.func.args(node.fun, "node.fun",
                                    function(x, labs, digits, varlen) NA)
        labs <- node.fun(x, labs, abs(digits), varlen)
        rpart.plot:::check.returned.labs(x, labs, node.fun.name)
    }
    labs <- paste0(prefix, labs, suffix)
    is.leaf <- rpart.plot:::is.leaf(frame)
    if(type == TYPE0.default || type == TYPE3.fancy.no.all)
        labs[!is.leaf] <- NA # no labels for internal nodes
    else if(type == TYPE5.varname.in.node) { # use split variable in interior nodes
        splits <- as.character(frame$var[!is.leaf])
        splits <- rpart.plot:::my.abbreviate(splits, varlen)
        labs[!is.leaf] <- splits
    }
    labs
}

#' This function is actually important and I actually modified it.
#' This might be what I want to make more customizable!!!
#' @keywords internal
#' @noRd
get.anova.labs <- function(x, extra, under, digits, xsep, varlen, under.percent)
{
    frame <- x$frame
    fitted <- frame$CI
    newline <- if(under) "\n\n" else "\n"
    ex <- if(extra < 100) extra else extra - 100
    labs <-
        if(ex == EX0)
            rpart.plot:::sprint("Fitted Mean:%s%s95%% CI:%s", rpart.plot:::format0(frame$y, digits), newline, fitted)
    else if(ex == EX1.NOBS) # add n?
        rpart.plot:::sprint("%s%sn=%s", fitted, newline, rpart.plot:::format0(frame$n, digits))
    else if (ex == EX2.CLASS.RATE) {
        #extra.help()
        stop("extra=", extra,
              ' is legal only for "class", "poisson" and "exp" models (you have an "anova" model)')
    }
    else if (ex > EX11.PROB.ACROSS.ALL.2ND.CLASS.DONT) {
        #extra.help()
        stop("extra=", extra, " is illegal")
    } else { # ex >= EX3.MISCLASS.RATE && ex <= EX11.PROB.ACROSS.ALL.2ND.CLASS.DONT
        #extra.help()
        stop("extra=", extra,
              ' is legal only for "class" models (you have an "anova" model)')
    }

    if(extra >= 100) {   # add percent?
        sep <- # space or newline before percent
            if(under.percent == 0) {
                "  "
            } else if(under.percent == 1) {
                if(under)
                    "\n\n"
                else
                    "\n"
            } else if(under.percent == 2) {
                if(extra == 100)
                    newline
                else
                    "  "
            }
        labs <- rpart.plot:::sprint("%s%s%s", labs, sep,
                       paste("n =", frame$wt))
    }
    labs
}
