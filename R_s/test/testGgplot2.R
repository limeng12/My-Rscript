library(ggplot2)
x <- seq(0, 6, .5)
y.a <- .1 * x -.1
y.b <- sin(x)
df <- data.frame(x=x, y=y.a, case='a')
df <- rbind(df, data.frame(x=x, y=y.b, case='b'))
print(ggplot(df) + geom_point(aes(x, y), color=ifelse(df$y<0, 'red', 'black')))