#!/bin/Rscript

######################################################################################
#
# Mike Mansfield 2021
# The input to this script is a Newick-format tree.
# The output is a PDF, optionally manually rooted (default)
# is midpoint-rooted.
#
######################################################################################

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(phangorn))
suppressPackageStartupMessages(require(phylotate))

####################################################################
#
# OPTPARSE AND INPUT OPTIONS
#
####################################################################
option_list = list(
	# Input options
	make_option(c('-t', '--treeFile'), action='store', default=NA, type='character',
		help='Path to a Newick-format tree file. Default: none'),
	make_option(c('-T', '--treeFormat'), action='store', default='newick', type='character',
		help='File format of input tree file. Default: newick'),
	make_option(c('-r', '--root'), action='store', default=NA, type='character',
		help='Tip to root tree by. Default: unrooted; can select a tip label or midpoint'),
	make_option(c('-s', '--sspb'), action='store', default=NA, type='character',
		help='Site-specific placement bias test result. Default: none, not plotted'),
	# Output options
	make_option(c('-o', '--outFile'), action='store', default=NA, type='character',
		help='Base name to use for output files. Default: [tree file base name]')
)

opt = parse_args(OptionParser(option_list=option_list))

# Make sure the tree file exists. If not, exit.
if(is.na(opt$treeFile)){
	cat('Error: an input tree must be specified with -t or --treeFile!\n', file=stderr())
	quit(status=1)
}

############################################################
# Output file checking.
############################################################
# If --outfile is not specified, use the tree file name.
if(is.na(opt$outFile)){
	output_base = basename(opt$treeFile)
} else{
	output_base = opt$outFile
}

############################################################
# Tree rooting
############################################################
# If rooting is not specified, use an unrooted tree.
if(is.na(opt$root)){
	root_type = 'unrooted'
} else {
	root_type = opt$root
}


############################################################
# Site-specific placement bias work
############################################################
if(is.na(opt$sspb)){
	sspb = FALSE
} else {
	sspb = TRUE
}

read_sspb <- function(sspb){
	df = read.delim(file=sspb, header=F, sep=' ')
	colnames(df) <- c('Position', 'Bias')
	df$Bias[which(df$Bias <0)] <- 0
	return(df)
}


############################################################
# Tree-reading
############################################################
read_tree_multiformat <- function(tree_file, tree_file_format='newick'){
	if( tree_file_format == "newick" ){
		tree = read.tree(tree_file)
	} else if( tree_file_format == "nexus" ){
		# Note that the phylotate's package is used here to parse MrBayes
		# better than the default ape read.nexus function.
		# I add the node labels manually.
		tree = read_annotated(tree_file, format='nexus')
		tree$node.labels = mb_attrs(tree)$prob*100
	}

	return(tree)
}


####################################################################
#
# CALCULATION AND PLOTTING FUNCTIONS
#
####################################################################
############################################################
# A helper function to make sure the root options make
# sense (the root option must be either a tip label or
# 'midpoint'.
############################################################
check_root_option <- function(root_option, tree_object){
	if(root_option %in% c('unrooted', 'midpoint', tree_object$tip.label)){
		return(TRUE)
	} else {
		return(FALSE)
	}
}

plot_sspb = function(sspb_df, sspb_base_name, col="brown3"){
	plot(sspb_df, main=paste(sspb_base_name, '\nSite-specific placement bias test result', sep=''), col=col, type='l' )
}

plot_tree = function(tree_object, root_option, show_node_label=T, tip_color='#595959', tree_base_name, align_tip_label=T){
	if(root_option=='unrooted') {
		plot(tree_object, type='unrooted', cex=0.66, tip.color=tip_color, show.node.label=show_node_label, align.tip.label=align_tip_label, main=paste(tree_base_name, '\nroot: unrooted', sep=''))
	} else if(root_option=='midpoint') {
		plot(midpoint(tree_object), cex=0.66, tip.color=tip_color, show.node.label=show_node_label, align.tip.label=align_tip_label, main=paste(tree_base_name, '\nroot: midpoint', sep=''))
	} else {
		plot(root(tree_object, outgroup=root_option), cex=0.66, tip.color=tip_color, show.node.label=show_node_label, align.tip.label=align_tip_label, main=paste(tree_base_name, '\nroot: ', root_option, sep=''))
	}
	add.scale.bar()
}


####################################################################
#
# WORKFLOW
#
####################################################################
if( sspb == TRUE ) {
	sspb_df = read_sspb(opt$sspb)
	pdf(paste(output_base, '.sspb.pdf', sep=''))
	plot_sspb(sspb_df, output_base, col="brown3")
	invisible(dev.off())
}

input_tree <- read_tree_multiformat(tree_file = opt$treeFile, tree_file_format=opt$treeFormat)

if(check_root_option(root_option=root_type, tree_object=input_tree)){
	pdf(paste(output_base, '.pdf', sep=''))
	plot_tree(tree_object=input_tree, root_option=root_type, show_node_label=T, tip_color='#595959', tree_base_name=output_base)
	invisible(dev.off())
} else {
	cat(paste('Error: the root ', root_type, ' is not valid!\nAvailable options: \'unrooted\', \'midpoint\', or any of the following:\n', paste(input_tree$tip.label, collapse=', ')), file=stderr())
	quit(status=1)
}
