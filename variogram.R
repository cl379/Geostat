data <- read.table("exp.csv", sep=";", header = TRUE, dec=".")

# semivariogram function. Needs a data.frame with $X, $Y and $var 
semiVariog <- function(data) {
    # cumulatives of pairing loop
    pitag <- c()
    sqdif <- c()
    # pairing loop
    for (i in 1:nrow(data)) {
	    for (j in 1:nrow(data)) {
		    if (j > i) {
                # pythagoras' theorem
			    pitag <- append(pitag, sqrt(((data[i, ]$X - data[j, ]$X)^2) + ((data[i, ]$Y - data[j, ]$Y)^2)))
                # square of differences
			    sqdif <- append(sqdif, (data[i, ]$var - data[j, ]$var)^2)
		    }
    	}
    }
    # build data.frame and order by distance
    pairs <- data.frame(pitag = pitag, sqdif = sqdif)
    pairs <- pairs[order(pairs$pitag), ]
    # calculate the relevant distance and subsample the pairs
    relevDist <- sqrt((max(data$X) - min(data$X))^2 + (max(data$Y) - min(data$Y))^2) / 3
    pairs <- pairs[ which(pairs$pitag < relevDist), ]
    # split classes by distances
    distGroup <- split(pairs, pairs$pitag)
    # cumulatives loops through breakes
    np <- c()
    dists <- c()
    gamma <- c()
    # loops through breaks
    for (i in distGroup) {
        # calculate the semivariance and establish number of pairs and distances tags 
	    gamma <- append(gamma, (1 / (2 * nrow(i))) * sum(i$sqdif))
	    np <- append(np, nrow(i))
	    dists <- append(dists, i$pitag[1])
    }
    # return the data.frame
    return(data.frame(np=np, dists=dists, gamma=gamma))
}
