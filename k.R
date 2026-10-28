# ---------------------------------------------------------------------------- #
# Diskrečių struktūrų statistinė analizė | Kontrolinis 1 (2025-11-06)           #
# Sprendimas. Reikia tik bazinio R + MASS (MASS ateina kartu su R).            #
# ---------------------------------------------------------------------------- #

rm(list = ls())
library(MASS)          # loglm (4.8 užduotis)


# ============================================================================ #
# UŽDAVINYS 1. 2x2 lentelės ryšio stiprumas                                    #
# ============================================================================ #

t1 <- matrix(c(13, 7,
                9, 11), nrow = 2, byrow = TRUE,
             dimnames = list(Grupe = c("A", "B"), Atsakymas = c("Taip", "Ne")))
addmargins(t1)

a <- t1[1, 1]; b <- t1[1, 2]; c <- t1[2, 1]; d <- t1[2, 2]
n <- sum(t1)

# Pirmiausia nepriklausomumo hipotezė (be Yates pataisos, kad phi = sqrt(chi2/n)).
chi1 <- chisq.test(t1, correct = FALSE)
chi1
fisher.test(t1)

# Phi koeficientas
phi <- (a * d - b * c) / sqrt((a + b) * (c + d) * (a + c) * (b + d))
phi
sqrt(chi1$statistic / n)                    # tas pats per chi^2

# Cramerio V (2x2 atveju V = |phi|)
V1 <- sqrt(chi1$statistic / (n * (min(dim(t1)) - 1)))
V1

# Yule Q
Q <- (a * d - b * c) / (a * d + b * c)
Q

# Šansų santykis
OR <- (a * d) / (b * c)
OR

# Išvada: phi ~ 0.20, Q ~ 0.39 -> silpnas ryšys (skalė iš paskaitos:
# [0;0.25] nėra, [0.25;0.5] silpnas). Chi^2 p ~ 0.21 > 0.05, t.y. ryšys
# statistiškai nereikšmingas -- su n = 40 imtimi to ir reikia tikėtis.


# ============================================================================ #
# UŽDAVINYS 2. Ranginių kintamųjų ryšys (Spirmeno r, Kendalo tau)              #
# ============================================================================ #

mat <- c(56, 75, 45, 71, 62, 64, 58, 80, 76, 61)
lt  <- c(66, 70, 40, 60, 65, 56, 59, 77, 67, 63)

rbind(Mat = rank(mat), LT = rank(lt))      # rangai (surištų rangų nėra)

# Spirmeno koeficientas pagal formulę r = 1 - 6*sum(d^2) / (n(n^2 - 1))
d2 <- (rank(mat) - rank(lt))^2
n2 <- length(mat)
r_s <- 1 - 6 * sum(d2) / (n2 * (n2^2 - 1))
r_s

cor(mat, lt, method = "spearman")           # patikrinimas

# Reikšmingumas: T = r*sqrt((n-2)/(1-r^2)) ~ t(n-2)
T2 <- r_s * sqrt((n2 - 2) / (1 - r_s^2))
p2 <- 2 * pt(-abs(T2), df = n2 - 2)
c(T = T2, p = p2)

cor.test(mat, lt, method = "spearman")      # tikslus p (mažai imčiai)
cor.test(mat, lt, method = "kendall")       # papildomai Kendalo tau

# Išvada: r ~ 0.67, p ~ 0.033 (t aproks.) / 0.039 (tikslus) < 0.05 -> vidutiniškai stiprus teigiamas
# monotoninis ryšys, statistiškai reikšmingas.


# ============================================================================ #
# UŽDAVINYS 3. Kendalo tau ir Z testas                                         #
# ============================================================================ #

g1 <- c(2, 5, 3, 1, 4, 6)
g2 <- c(1, 6, 4, 3, 2, 5)
n3 <- length(g1)

# Suderintų (nc) ir nesuderintų (nd) porų skaičius
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
cor(g1, g2, method = "kendall")             # patikrinimas

# Z = 3*tau*sqrt(n(n-1)) / sqrt(2(2n+5)) ~ N(0,1)
# (lapelyje vardiklyje šaknis nematoma, bet teisinga formulė su šaknimi --
#  ji ekvivalenti paskaitos Z = tau / sqrt(2(2n+5) / (9n(n-1))).)
Z3 <- 3 * tau * sqrt(n3 * (n3 - 1)) / sqrt(2 * (2 * n3 + 5))
p3 <- 2 * pnorm(-abs(Z3))
c(Z = Z3, z_kr = qnorm(0.975), p = p3)

cor.test(g1, g2, method = "kendall")        # tikslus testas palyginimui

# Išvada: nc = 11, nd = 4, tau = 0.467, |Z| = 1.32 < 1.96,
# p ~ 0.19 (tikslus p = 0.27) -> H0: tau = 0 neatmetama. Gydytojų nuomonės iš dalies sutampa,
# bet su 6 pacientais negalima teigti, kad jos statistiškai reikšmingai suderintos.


# ============================================================================ #
# UŽDAVINYS 4. Titanic (R įdiegta 4-matė lentelė: Class, Sex, Age, Survived)   #
# ============================================================================ #

str(Titanic)
dimnames(Titanic)

# --- 4.1 Class x Survived su marginaliniais dažniais ------------------------- #
cs <- margin.table(Titanic, c("Class", "Survived"))
addmargins(cs)
round(prop.table(cs, margin = "Class"), 3)   # išgyvenamumas pagal klasę

# Išvada: iš 2201 žmonių išgyveno 711 (32.3 %). Daugiausia žuvo įgulos (673)
# ir 3 klasės (528). Išgyvenusiųjų dalis krenta nuo 1 klasės (62 %) iki 3
# klasės (25 %) ir įgulos (24 %).

# --- 4.2 Lyties dažnių lentelė ---------------------------------------------- #
sex <- margin.table(Titanic, "Sex")
sex
names(which.max(sex))
# Vyrų 1731, moterų 470 -> vyrų beveik 4 kartus daugiau.

# --- 4.3 Sex x Survived ir sąlyginės tikimybės ------------------------------- #
ss <- margin.table(Titanic, c("Sex", "Survived"))
addmargins(ss)
p_cond <- prop.table(ss, margin = "Sex")     # P(Survived | Sex)
round(p_cond, 3)
# P(išgyveno | vyras) = 0.212, P(išgyveno | moteris) = 0.732

# --- 4.4 Šansai ir šansų santykis ------------------------------------------- #
odds_m <- ss["Male", "Yes"] / ss["Male", "No"]
odds_f <- ss["Female", "Yes"] / ss["Female", "No"]
c(sansas_vyru = odds_m, sansas_moteru = odds_f)

OR_mf <- odds_m / odds_f
OR_fm <- odds_f / odds_m
c(OR_vyrai_vs_moterys = OR_mf, OR_moterys_vs_vyrai = OR_fm)
# Vyrui šansas išgyventi 0.27 (t.y. ~1 išgyvenęs 3.7 žuvusiems).
# OR ~ 0.10: vyrų šansai išgyventi ~10 kartų mažesni nei moterų.

# --- 4.5 Nepriklausomumo testas Sex x Survived ------------------------------ #
chisq.test(ss)                       # su Yates pataisa (2x2)
chisq.test(ss, correct = FALSE)
fisher.test(ss)                      # duoda ir OR su 95 % PI
summary(ss)                          # kaip paskaitoje (summary(table))
# p < 2.2e-16 -> H0 apie lyties ir išgyvenamumo nepriklausomumą atmetama.

# --- 4.6 Grafikas Class x Survived ------------------------------------------ #
op <- par(mfrow = c(1, 2))
mosaicplot(cs, color = c("tomato", "seagreen3"),
           main = "Class x Survived (mozaika)")
barplot(t(prop.table(cs, margin = 1)), beside = FALSE,
        col = c("tomato", "seagreen3"), legend.text = TRUE,
        args.legend = list(x = "topright", bty = "n", cex = 0.8),
        ylab = "Dalis", main = "Išgyvenusių dalis pagal klasę")
par(op)
# Mozaikoje stulpelio plotis = klasės dydis, žalios dalies aukštis =
# išgyvenusiųjų dalis. Matosi: įgula didžiausia grupė, o išgyvenamumas
# mažėja nuo 1 klasės prie 3 klasės ir įgulos.

# --- 4.7 Sex x Age x Survived; moterų dalis tarp (ne)išgyvenusių ------------ #
sas <- margin.table(Titanic, c("Sex", "Age", "Survived"))
ftable(sas, row.vars = c("Age", "Survived"))

# Moterų santykinis dažnis kiekvienoje (Age, Survived) grupėje:
p_sex <- prop.table(sas, margin = c("Age", "Survived"))
round(ftable(p_sex, row.vars = c("Age", "Survived")), 3)
round(p_sex["Female", , ], 3)
# Suaugę: tarp išgyvenusių moterų 0.483, tarp žuvusių 0.076.
# Vaikai: tarp išgyvenusių moterų 0.491, tarp žuvusių 0.327.

# --- 4.8 Ar Survived ir Age sąlyginai nepriklausomi, kai žinoma Sex? -------- #
# H0: Age ⟂ Survived | Sex.

# (a) Atskiri chi^2 testai kiekviename lyties sluoksnyje
for (s in dimnames(sas)$Sex) {
  cat("\n---", s, "---\n")
  print(sas[s, , ])
  print(chisq.test(sas[s, , ]))
}

# (b) Cochran-Mantel-Haenszel testas (2x2xK lentelės, sluoksnis = Sex)
sas_cmh <- aperm(sas, c("Age", "Survived", "Sex"))
mantelhaen.test(sas_cmh)

# (c) Logtiesinis modelis [Sex*Age][Sex*Survived] vs prisotintas (G^2 testas,
#     kaip 11 temoje). Šis modelis kaip tik ir yra sąlyginio nepriklausomumo.
m_ci  <- loglm(~ Sex * Age + Sex * Survived, data = sas)
m_ci
m_sat <- loglm(~ Sex * Age * Survived, data = sas)
anova(m_ci, m_sat)

# Išvada: CMH p = 0.019, logtiesinio modelio G^2 = 22.2, df = 2, p = 1.5e-5
# -> H0 atmetama. Net žinant lytį, amžius ir
# išgyvenamumas priklausomi -- vaikų išgyvenamumas didesnis (ypač berniukų
# palyginti su suaugusiais vyrais). Pastaba: moterų sluoksnyje ryšys silpnas
# ir nereikšmingas, visą efektą neša vyrų sluoksnis.