x <- 1:26
xboot <- sample(x, replace = TRUE)
xboot

y <- LETTERS
y[xboot]
