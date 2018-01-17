import interaction.util as iu
import interaction.plane_regression as plane


def calc_pi_pi(ring1, ring2, aromDist=6, aromDihedralAng=30, aromDispAng=20):

    distance = iu.calc_euclidean_distance(ring1.centroid, ring2.centroid)

    data = []
    if (distance <= aromDist):

        ring1.calc_normal()
        ring2.calc_normal()

        dihedralAngle = iu.to_quad1(iu.calc_angle(ring1.normal, ring2.normal))

        vectorCC = ring2.centroid - ring1.centroid
        dispAngle1 = iu.to_quad1(iu.calc_angle(ring1.normal, vectorCC))
        dispAngle2 = iu.to_quad1(iu.calc_angle(ring2.normal, vectorCC))

        normal1 = plane.get_orthog_point(ring1.centroid, ring1.normal, 0.05)
        normal2 = plane.get_orthog_point(ring2.centroid, ring2.normal, 0.05)

        interaction = [0, 0, 0]
        if (dihedralAngle > aromDihedralAng):
            interaction = [0, 0, 1]
        else:
            if (dispAngle1 <= aromDispAng):
                interaction = [1, 0, 0]
            else:
                interaction = [0, 1, 0]

        data += interaction \
            + [distance, dihedralAngle, dispAngle1, dispAngle2] \
            + list(ring1.centroid) + list(normal1) \
            + list(ring2.centroid) + list(normal2)

    return data
