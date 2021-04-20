"""
Предварительный расчёт ЦБК
Проверка основных параметров
"""

import numpy as np

pi = np.pi
print(pi)
grad_rad = pi / 180


def eps_gdf(lam, kk=1.4):
    return (1 - (kk - 1) / (kk + 1) * lam ** 2) ** (1 / (kk - 1))


def a_kr(tt_0_t, rr=287, kk=1.4):
    return np.sqrt(2 * kk / (kk - 1) * rr * tt_0_t)


# Ввод параметров

G = 4
pika_t = 4
p_0_t = 101325
t_0_t = 287
beta_2l = 60

# Свойства рабочего тела
R = 287
k = 1.4
Cp = 1004


# Подбор параметров
z_opt = 24
D_ = 0.55    # относительный диаметр периферии при входе D1п / D2 (0.45 - 0.65)
d_ = 0.45    # относительный диамтр втулки при входе (0.35 - 0.55)
C_2m0_ = 0.31  # обобшённый коэффициент расхода (0.22 - 0.4)
eta_p_t_k = 0.84
alfa_1_sr = 90  # угол закрутки лопатки при входе
gamma_1_sr = 5  # угол наклона средней линии тока при входе
k_cm = 1  # коэффициент разгона меридиональной скорости
sigma_vy = 0.99  # коэффициент сохранения полного давления во входном патрубке
#  sigma_vy = (0.985 - 0.995)
mu_vy = 0.99  # коэф. загромождения при входе (0.98 - 0.99)

H_k0_ = 0.9  # коэф. эффективного напора (0.87 - 0.93)
eta_ad_k = (pika_t ** ((k - 1) / k) - 1) / (pika_t ** ((k - 1) / (k * eta_p_t_k))-1)

H_k_ = H_k0_ * np.sqrt(np.sin(beta_2l * grad_rad))

H_ad_t = Cp * t_0_t * (pika_t ** ((k - 1) / k) - 1)
H_ = H_k_ * eta_ad_k
U_2 = np.sqrt(H_ad_t / H_)
U_1p = U_2 * D_
C_2m_ = C_2m0_ * np.sin(beta_2l * grad_rad)
C_2m = C_2m_ * U_2
C_1m = C_2m * k_cm
C_1 = C_1m / np.sin(alfa_1_sr * grad_rad)
lam_1 = C_1 / a_kr(t_0_t, R, k)

ro_1_t = p_0_t * sigma_vy / (R * t_0_t)

F_1 = G / (C_1m * ro_1_t * eps_gdf(lam_1) * mu_vy)

print('eta_ad_k = ', eta_ad_k)
print('U_2 = ', U_2)
# Расчёт
