rm(list = ls())

# 1

t1 <- matrix(c(13, 7,
                9, 11), nrow = 2, byrow = TRUE,
             dimnames = list(Grupe = c("A", "B"), Atsakymas = c("Taip", "Ne")))
addmargins(t1)

a <- t1[1, 1]; b <- t1[1, 2]; c <- t1[2, 1]; d <- t1[2, 2]
n <- sum(t1)

chi1 <- chisq.test(t1, correct = FALSE)
chi1
fisher.test(t1)

phi <- (a * d - b * c) / sqrt((a + b) * (c + d) * (a + c) * (b + d))
phi
sqrt(chi1$statistic / n)

V1 <- sqrt(chi1$statistic / (n * (min(dim(t1)) - 1)))
V1

Q <- (a * d - b * c) / (a * d + b * c)
Q

OR <- (a * d) / (b * c)
OR

# Rezultatai: phi = V = 0.20, Q = 0.39, OR = 2.27, chi^2 = 1.62, p = 0.20.
# Išvada: silpnas ryšys, statistiškai nereikšmingas (p > 0.05).


# 2

mat <- c(56, 75, 45, 71, 62, 64, 58, 80, 76, 61)
lt  <- c(66, 70, 40, 60, 65, 56, 59, 77, 67, 63)

rbind(Mat = rank(mat), LT = rank(lt))

d2  <- (rank(mat) - rank(lt))^2
n2  <- length(mat)
r_s <- 1 - 6 * sum(d2) / (n2 * (n2^2 - 1))
r_s
cor(mat, lt, method = "spearman")

T2 <- r_s * sqrt((n2 - 2) / (1 - r_s^2))
p2 <- 2 * pt(-abs(T2), df = n2 - 2)
c(T = T2, kritine = qt(0.975, n2 - 2), p = p2)

cor.test(mat, lt, method = "spearman")

# Rezultatai: r = 0.673, T = 2.57 > t(8) = 2.306, p = 0.033.
# Išvada: vidutiniškai stiprus teigiamas monotoninis ryšys, statistiškai
# reikšmingas.


# 3

g1 <- c(2, 5, 3, 1, 4, 6)
g2 <- c(1, 6, 4, 3, 2, 5)
n3 <- length(g1)

nc <- 0; nd <- 0
for (i in 1:(n3 - 1)) {
  for (j in (i + 1):n3) {
    s <- sign(g1[i] - g1[j]) * sign(g2[i] - g2[j])
    if (s > 0) nc <- nc + 1
    if (s < 0) nd <- nd + 1
  }
}
c(nc = nc, nd = nd)

tau <- (nc - nd) / (n3 * (n3 - 1) / 2)
tau
cor(g1, g2, method = "kendall")

Z3 <- 3 * tau * sqrt(n3 * (n3 - 1)) / sqrt(2 * (2 * n3 + 5))
p3 <- 2 * pnorm(-abs(Z3))
c(Z = Z3, kritine = qnorm(0.975), p = p3)

# Rezultatai: nc = 11, nd = 4, tau = 0.467, |Z| = 1.32 < 1.96, p = 0.19.
# Išvada: H0 neatmetama. Gydytojų nuomonės iš dalies sutampa, bet
# statistiškai reikšmingai suderintomis jų laikyti negalima.


# 4.1

cs <- margin.table(Titanic, c("Class", "Survived"))
addmargins(cs)
round(prop.table(cs, margin = "Class"), 3)

# Rezultatai: iš 2201 išgyveno 711 (32.3 %). Išgyvenusių dalis: 1 klasė 62.5 %,
# 2 klasė 41.4 %, 3 klasė 25.2 %, įgula 24.0 %.
# Išvada: kuo žemesnė klasė, tuo mažesnis išgyvenamumas; daugiausia žuvo
# įgulos (673) ir 3 klasės (528).


# 4.2

sex <- margin.table(Titanic, "Sex")
sex
names(which.max(sex))

# Rezultatai: vyrų 1731, moterų 470.
# Išvada: vyrų beveik 4 kartus daugiau.


# 4.3

ss <- margin.table(Titanic, c("Sex", "Survived"))
addmargins(ss)
round(prop.table(ss, margin = "Sex"), 3)

# Rezultatai: P(išgyveno | vyras) = 0.212, P(išgyveno | moteris) = 0.732.


# 4.4

odds_m <- ss["Male", "Yes"] / ss["Male", "No"]
odds_f <- ss["Female", "Yes"] / ss["Female", "No"]
c(sansas_vyru = odds_m, sansas_moteru = odds_f)

c(OR_vyrai_vs_moterys = odds_m / odds_f, OR_moterys_vs_vyrai = odds_f / odds_m)

# Rezultatai: vyrų šansas 0.27, moterų 2.73, OR = 0.099.
# Išvada: vyrų šansai išgyventi ~10 kartų mažesni nei moterų.


# 4.5

O <- ss
n <- sum(O)
E <- outer(rowSums(O), colSums(O)) / n
E

X2    <- sum((O - E)^2 / E)
df    <- (nrow(O) - 1) * (ncol(O) - 1)
kr    <- qchisq(0.95, df)
p_val <- pchisq(X2, df, lower.tail = FALSE)
c(X2 = X2, df = df, kritine = kr, p = p_val)

# Rezultatai: chi^2 = 456.9 > 3.84, df = 1, p ~ 0.
# Išvada: H0 atmetama, lytis ir išgyvenamumas priklausomi.


# 4.6

barplot(t(cs), beside = TRUE,
        col = c("tomato", "seagreen3"),
        legend.text = c("Žuvo", "Išgyveno"),
        args.legend = list(x = "topleft", bty = "n"),
        xlab = "Klasė", ylab = "Žmonių skaičius",
        main = "Išgyvenusieji ir žuvusieji pagal klasę")

# Išvada: tik 1 klasėje išgyvenusių daugiau nei žuvusių; 3 klasėje ir įguloje
# žuvusių kelis kartus daugiau.


# 4.7

sas <- margin.table(Titanic, c("Sex", "Age", "Survived"))
ftable(sas, row.vars = c("Age", "Survived"))

p_sex <- prop.table(sas, margin = c("Age", "Survived"))
round(ftable(p_sex, row.vars = c("Age", "Survived")), 3)
round(p_sex["Female", , ], 3)

# Rezultatai: suaugę: moterų dalis tarp išgyvenusių 0.483, tarp žuvusių 0.076.
# Vaikai: tarp išgyvenusių 0.491, tarp žuvusių 0.327.
# Išvada: tarp išgyvenusių moterų dalis daug didesnė, ypač suaugusių grupėje.


# 4.8

X2_k <- c()
for (s in dimnames(sas)$Sex) {
  O_k <- sas[s, , ]
  E_k <- outer(rowSums(O_k), colSums(O_k)) / sum(O_k)
  X2_k[s] <- sum((O_k - E_k)^2 / E_k)
  cat("\n---", s, "---\n")
  print(O_k)
  cat("chi^2 =", round(X2_k[s], 3),
      " p =", signif(pchisq(X2_k[s], 1, lower.tail = FALSE), 3), "\n")
}

X2_bendra <- sum(X2_k)
df8       <- length(X2_k) * 1
c(X2 = X2_bendra, df = df8,
  kritine = qchisq(0.95, df8),
  p = pchisq(X2_bendra, df8, lower.tail = FALSE))

# Rezultatai: vyrai chi^2 = 23.1 (p = 1.5e-6), moterys chi^2 = 3.05 (p = 0.08).
# Bendra chi^2 = 26.2 > 5.99, df = 2, p = 2e-6.
# Išvada: H0 atmetama. Net žinant lytį, amžius ir išgyvenamumas priklausomi.
# Efektą neša vyrai (berniukai išgyveno daug dažniau nei suaugę vyrai);
# moterų grupėje ryšys nereikšmingas.
