import pygame
import math

two_pi = 2 * math.pi
gravity = None


class homogeneousRod:
    def __init__(self, cm, l, d, theta, mass=1):
        self.mass = mass
        self.inv_mass = 1/mass
        self._cm = cm
        self._fixed_point = (0, 0)
        self._targeted = False
        self.cmV = (0, 0)
        self.cmA = (0, 0)
        self._l = l
        self._d = d
        if d != 0 and l != 0:
            self.inv_inertMom = 12 / ( mass * ((2 * l)**2 + (2 * d)**2) )
        else: self.inv_inertMom = 0
        self._theta = theta
        self.thetaV = 0
        self.thetaA = 0
        self._rv = (math.cos(theta), math.sin(theta))
        self._tv = (self._rv[1], -self._rv[0])
        self.update_vertexes()

    @property
    def cm(self):
        return self._cm

    @property
    def vertexes(self):
        return self._vertexes

    @property
    def fixed_point(self):
        return self._fixed_point

    @property
    def targeted(self):
        return self._targeted

    def target(self):
        self._targeted = not self._targeted

    def update_vertexes(self):
        l1 = (self._l * self._rv[0], self._l * self._rv[1])
        l2 = (- self._l * self._rv[0], - self._l * self._rv[1])
        d1 = (self._d * self._tv[0], self._d * self._tv[1])
        d2 = (- self._d * self._tv[0], - self._d * self._tv[1])
        self._vertexes = [(self._cm[0] + l1[0] + d1[0], self._cm[1] + l1[1] + d1[1]),
                          (self._cm[0] + l2[0] + d1[0], self._cm[1] + l2[1] + d1[1]),
                          (self._cm[0] + l2[0] + d2[0], self._cm[1] + l2[1] + d2[1]),
                          (self._cm[0] + l1[0] + d2[0], self._cm[1] + l1[1] + d2[1])]

    @property
    def tv(self):
        return self._tv

    @property
    def rv(self):
        return self._rv

    def update_attributes(self, l=None, d=None, mass=None):
        if mass is not None:
            self.mass = mass
            self.inv_mass = 1/mass
        if l is not None:
            self._l = l
        if d is not None:
            self._d = d

        if ( mass is not None or l is not None or d is not None ) and l + d != 0:
            self.inv_inertMom = 12 / ( self.mass * ((2 * self._l)**2 + (2 * self._d)**2) )

        if l is not None or d is not None:
            self.update_vertexes()

    def fix_point(self, point):
        rel_point = (point[0] - self._cm[0], point[1] - self._cm[1])
        self._fixed_point = (rel_point[0] * self.rv[0] + rel_point[1] * self.rv[1],
                             rel_point[0] * self.tv[0] + rel_point[1] * self.tv[1])

    def update_angle(self, newtheta):
        self._theta = newtheta
        self._rv = (math.cos(self._theta), math.sin(self._theta))
        self._tv = (self._rv[1], -self._rv[0])
        self.update_vertexes()

    def update_cm(self, newcm):
        self._cm = newcm
        self.update_vertexes()

    def update_kinematics(self, dt):
        self.thetaV += self.thetaA * dt
        if self.thetaV != 0:
            if ( na := self._theta + dt * self.thetaV ) >= two_pi: na += - two_pi
            elif na < 0: na += two_pi
            self._theta = na
            self._rv = (math.cos(self._theta), math.sin(self._theta))
            self._tv = (self._rv[1], -self._rv[0])

        self.cmV = (self.cmV[0] + dt * self.cmA[0], self.cmV[1] + dt * self.cmA[1])
        self._cm = (self.cm[0] + dt * self.cmV[0], self.cm[1] + dt * self.cmV[1])
        self.update_vertexes()

        self.thetaA = 0
        self.cmA = (0, 0)

    def collision_point(self, point):
        relpos = (point[0] - self._cm[0], point[1] - self._cm[1])
        t_component = abs(relpos[0] * self._tv[0] + relpos[1] * self._tv[1])
        r_component = abs(relpos[0] * self._rv[0] + relpos[1] * self._rv[1])
        if r_component <= self._l and t_component <= self._d:
            return True
        else: return False

    def handle_forces(self, F, contact_point):
        vcp = (contact_point[0] - self._cm[0], contact_point[1] - self._cm[1])
        self.thetaA += (vcp[0] * F[1] - vcp[1] * F[0]) * self.inv_inertMom
        self.cmA = (self.cmA[0] + F[0] * self.inv_mass, self.cmA[1] + F[1] * self.inv_mass)

def force_calculator(vx, vy, dt, vcm_x, vcm_y, r_cm_p_x, r_cm_p_y, ang_vel, inv_inert_mom, inv_mass):
    h = dt * inv_inert_mom * r_cm_p_x * r_cm_p_y
    k_y = dt * ( inv_mass + (r_cm_p_y ** 2) * inv_inert_mom )
    k_x = dt * ( inv_mass + (r_cm_p_x ** 2) * inv_inert_mom )
    if ( hkxy := h**2 - k_y * k_x ) == 0:
        print("h = k_y")
        return (0, 0)
    F_y = ( k_y * (-vy + vcm_y + ang_vel * r_cm_p_x) + h * (vcm_x - vx - ang_vel * r_cm_p_y) ) / hkxy
    F_x = ( k_x * (-vx + vcm_x - ang_vel * r_cm_p_y) + h * (-vy + vcm_y + ang_vel * r_cm_p_x) ) / hkxy

    return (F_x, F_y)

try:
    gravity = int(input('force of gravity:\n'))
except (ValueError):
    while True:
        try:
            n = int(input('type a number!!\n'))
            break
        except(ValueError):
            continue

pygame.init()
clock = pygame.time.Clock()
WIDTH = 1000
HEIGHT = 800
dt = 0.001
inv_dt = 1/dt
screen = pygame.display.set_mode((WIDTH, HEIGHT))

rods = []
run = True
buffer_rod = None
setting_rod = False
anglebuff = 0
last_mp = (0, 0)
actual_mp = (0, 0)
interaction_type = 0
fixated = False

while run:
    screen.fill((0, 0, 0))
    last_mp = actual_mp
    actual_mp = pygame.mouse.get_pos()

    for event in pygame.event.get():

        if event.type == pygame.MOUSEMOTION:
            mp = pygame.mouse.get_pos()
            if setting_rod:
                pv = (mp[0] - buffer_rod.cm[0], mp[1] - buffer_rod.cm[1])
                d = abs(buffer_rod.tv[0] * pv[0] + buffer_rod.tv[1] * pv[1])
                l = abs(buffer_rod.rv[0] * pv[0] + buffer_rod.rv[1] * pv[1])
                buffer_rod.update_attributes(l=l, d=d)

        if event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1:
                buffer_rod = homogeneousRod(pygame.mouse.get_pos(), 0, 0, anglebuff, mass=3)
                setting_rod = True

            elif event.button == 2:
                fixated = not fixated
                if fixated:
                    for rod in rods:
                        rod.fix_point(pygame.mouse.get_pos())
                        rod.target()
                else:
                    for rod in rods:
                        rod.target()

            elif event.button == 3:
                interaction_type = (interaction_type + 1) % 2

        elif event.type == pygame.MOUSEWHEEL:
            if ( na := event.y * 0.1 + anglebuff ) >= 2 * math.pi: na += - two_pi
            elif na < 0: na += two_pi

            anglebuff = na
            if buffer_rod is not None:
                buffer_rod.update_angle(anglebuff)

            elif rods != []:
                rods[-1].update_angle(anglebuff)

        elif event.type == pygame.MOUSEBUTTONUP:
            if setting_rod:
                setting_rod = False
                rods.append(buffer_rod)
                buffer_rod = None

        elif event.type == pygame.KEYDOWN:
            if event.key == pygame.K_d:
                if buffer_rod is not None:
                    buffer_rod = None
                    setting_rod = False

                elif rods != []:
                   del rods[-1]

        elif event.type == pygame.QUIT:
            run = False


    v = None
    for rod in rods:

        if interaction_type != 0 and ( rod.collision_point(actual_mp) or rod.targeted ):
            rel_pos = (actual_mp[0] - rod.cm[0], actual_mp[1] - rod.cm[1])
            if v == None: v = ( (actual_mp[0] - last_mp[0]) * inv_dt, (actual_mp[1] - last_mp[1]) * inv_dt)
            F = force_calculator(v[0], v[1], dt, rod.cmV[0], rod.cmV[1], rel_pos[0], rel_pos[1], rod.thetaV,
                                 rod.inv_inertMom, rod.inv_mass)
            rod.handle_forces(F, actual_mp)

        rod.handle_forces((0, gravity), rod.cm)
        rod.update_kinematics(dt)

        if rod.targeted:
            rv_c = rod.fixed_point[0]
            tv_c = rod.fixed_point[1]
            rod.update_cm((actual_mp[0] - rv_c * rod.rv[0] - tv_c * rod.tv[0], actual_mp[1] - rv_c * rod.rv[1] -
                            tv_c * rod.tv[1]))

        vert = rod.vertexes
        pygame.draw.line(screen, (255, 255, 255), vert[0], vert[1])
        pygame.draw.line(screen, (255, 255, 255), vert[1], vert[2])
        pygame.draw.line(screen, (255, 255, 255), vert[2], vert[3])
        pygame.draw.line(screen, (255, 255, 255), vert[3], vert[0])
        pygame.draw.circle(screen, (255, 255, 255), rod.cm, 1)

    if buffer_rod is not None:
        vert = buffer_rod.vertexes
        pygame.draw.line(screen, (255, 255, 255), vert[0], vert[1])
        pygame.draw.line(screen, (255, 255, 255), vert[1], vert[2])
        pygame.draw.line(screen, (255, 255, 255), vert[2], vert[3])
        pygame.draw.line(screen, (255, 255, 255), vert[3], vert[0])


    pygame.display.update()
    clock.tick(60)
