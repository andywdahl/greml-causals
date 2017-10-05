rm( list=ls() )
library(lme4)

klist		<- c( 2, 5, 10, 20, 50, 100 )
h2s		<- c( .25, .5, .75 )
phi0s	<- h2s / (1-h2s)

p			<- 1e5
ns		<- c( 2e3, 1e4 )
seeds	<- c( 696, 969 )

xs			<- 2:max(klist)
theor	<- function(k,wb,t,n,h2)
  2*(1-wb)^2*h2^2*(1-h2)^2/(t*k*n)

pdf( 'Figure3.pdf', width=10, height=7 )
par( mar=c(5,5,1,1) )

#plot(	c(1,max(klist)), c(0,sqrt(.003)), type='n', axes=F, xlab='# of Causal SNPs', ylab=expression(paste(h["*"]^2,' Standard Deviation')) )
plot(	c(1,max(klist)), c(0,.05), type='n', axes=F, xlab='# of Causal SNPs', ylab=expression(paste(h["*"]^2,' Standard Deviation')) )
box()
axis( 1, at=klist+c(-.2,.2,rep(0,length(klist)-2)), lab=klist )
axis(2)
abline( h=0, col='lightgrey', lwd=2, lty=1 )

leg=c( expression(paste( h^2, '=', 25, '%, n=2,000' )),
			expression(paste( h^2, '=', 50, '%, n=2,000' )),
			expression(paste( h^2, '=', 75, '%, n=2,000' )),
			expression(paste( h^2, '=', 25, '%, n=10,000' )),
			expression(paste( h^2, '=', 50, '%, n=10,000' )),
			expression(paste( h^2, '=', 75, '%, n=10,000' )) )
legend( 'topright', leg=leg,
	lty=c(1	,1	,1	,2,2,2),
	pch=c(16,16	,16	,1,1,1),
	col=c(1	,2	,3	,1,2,3) )

for( n.i in seq( along=ns ) )
	for( seed.i in seq( along=seeds ) )
{
	load( paste0( 'Rdata/Zfile_'		, ns[n.i], '_', p, '_', seeds[seed.i], '.Rdata' ) )## kinship eigenvals
	load( paste0( 'Rdata/savefile_'	, ns[n.i], '_', p, '_', seeds[seed.i], '.Rdata' ) )## estimated phis

	v		<- apply( allphi / ( 1 + allphi ), 1:2, function(h2){
		jj=rep(1:nrow(h2),ncol(h2))
		VarCorr( lmer(c(h2)~1|jj) )[[1]][1]
	})

	for( h.i in seq( along=h2s ) ){
		w						<- 1/(1+phi0s[h.i]*s2)
		v_theor			<- theor( xs, wb=mean(w), t=var(w), ns[n.i], h2s[h.i] )
		lines(	xs		, sqrt( v_theor	)	, col=h.i, lty=n.i )
		points( klist	, sqrt( v[h.i,])	, col=h.i, pch=ifelse( n.i==1, 16, 1 ), cex=1.2 )
	}
}
dev.off()
