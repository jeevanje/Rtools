Rsquared = function(x,y){
	# y = observed
	model = lm( as.vector(y) ~ as.vector(x) )
	return( round(summary(model)$r.squared,digits=2) )
} 