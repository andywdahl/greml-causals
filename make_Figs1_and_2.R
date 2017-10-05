rm( list=ls() )
set.seed( 5952 )

require(RMTstat)
n			<- 4000
p			<- 3e5
phi		<- 1
theta	<- 1

### MP spectrum
s2	<- rmp(n,svr=p/n)
z		<- rnorm(n,sd=sqrt(1+phi*s2))
Lz	<- log(z^2)


abfxn	<- function( phi, bsig=1, ... ){
	aa	<- log(1+phi*bsig)-phi*bsig/(1+bsig*phi)
	bb	<- phi/(1+bsig*phi)
	ltheta	<- -log(sum(z^2/(1+s2*phi))/n)
	abline(aa-ltheta,bb,...)
}

xlab	<- expression(s[i]^2)
ylab	<- expression(log(z[i]^2))

par( mar=c(8,7,1,1) )

pdf( 'Figure1.pdf' )
par( mar=c(5,5,5,5) )
plot(s2,Lz,ylab=ylab,xlab=xlab, main='')
abfxn( phi, col=2, lwd=3 )
for( tst in c( .5, 0, 2, 5 ) )
	abfxn( tst, col=3, lwd=3, lty=3 )
dev.off()

pdf( 'Figure2.pdf' )
par( mar=c(5,5,5,5) )
plot(s2,Lz,,ylab=ylab,xlab=xlab,ylim=0:1,main='')
abfxn( phi, col=2, lwd=3 )
for( tst in c( .5, 0, 2, 5 ) )
	abfxn( tst, col=3, lwd=3, lty=3 )
dev.off()
