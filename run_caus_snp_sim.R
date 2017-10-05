rm( list=ls() )

p				<- 1e5
klist		<- c( 2, 5, 10, 20, 50, 100 )
h2s			<- c( .25, .5, .75 )
numsig	<- length(h2s)
numgen	<- 1000
numphen	<- 100

phifxn	<- function(z,s2,lowlim=-.2,upplim=1e3){
	interior	<- score_eq(lowlim,z,s2)*score_eq(upplim,z,s2)<0
	if( ! interior )
		return( upplim )
	uniroot(score_eq,c(lowlim,upplim),z=z,s2=s2)$root
}
score_eq	<- function(phi,z,s2)
	-1/phi * cov( 1/(1+phi*s2), z^2/(1+phi*s2) )

for( n in c( 2e3, 1e4 ) )
	for( seed in sample(c( 696, 969 )) )
try({

	savefile	<- paste0( 'Rdata/savefile_', n, '_', p, '_', seed, '.Rdata' )
	if( file.exists( savefile ) ) next
	print( c( n, seed ) )

	## create/load Z + its svd
	Zfile			<- paste0( 'Rdata/Zfile_', n, '_', p, '_', seed, '.Rdata' )
	if( file.exists( Zfile ) ){
		load( Zfile )
	} else {
		set.seed(seed)
		mafs	<- runif(p,.05,.95)
		Z			<- scale( sapply( mafs, function(p0) rbinom(n,2,p0) ) )
		s			<- svd(Z/sqrt(p))
		s2		<- s$d^2
		UZ	<- t(s$u)%*%Z
		save( UZ, mafs, s, s2, file=Zfile )
	}

	## run many h2 simulations for fixed Z
	allphi	<- array( NA, dim=c( length(h2s), length(klist), numgen, numphen ) )
	for( k in klist )
		for( sig.i in seq( along=h2s ) )
			for( i in 1:numgen )
	{
		if( i == 1 ) print( c( k, sig.i ) )
		UZc	<- UZ[,sample(p,k),drop=F]
		for(j in 1:numphen){

			## draw u, fix size to give right h2
			u		<- rnorm(k)
			u		<- u/sqrt(sum(u^2))

			## create pheno
			z	<- sqrt(h2s[sig.i]) * UZc %*% u + sqrt(1-h2s[sig.i]) * rnorm(n)

			## run h2 on whitened phenos
			allphi[sig.i,which(klist==k),i,j]	<- phifxn( z=z, s2 )
		}
	}
	save( allphi, file=savefile )
})
