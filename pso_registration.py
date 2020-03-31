import numpy, pygmo, open3d, copy, argparse, os, logging, tqdm, csv, sys
from scipy.spatial.transform import Rotation
import scipy.spatial
from math import sin, cos, radians, exp
import scipy.stats
from numba import jit, float64


class dummy_isl:
    def run_evolve(self, algo, pop):
        new_pop = algo.evolve(pop)
        return algo, new_pop

    def get_name(self):
        return "Dummy slow island"


class pso_registration:
    def __init__(self, source, target, bounds):
        self._source = copy.deepcopy(source)
        self._target = copy.deepcopy(target)
        self._target_tree = scipy.spatial.cKDTree(self._target.points, leafsize=10)
        self._source_centroid, _ = self._source.compute_mean_and_covariance()
        self._bounds = bounds

    def spherical_to_rotvec(inclination, radius, azimuth):
        x = radius * sin(radians(inclination)) * cos(radians(azimuth))
        y = radius * sin(radians(inclination)) * sin(radians(azimuth))
        z = radius * cos(radians(inclination))
        vec = [x, y, z]
        return Rotation.from_rotvec(vec)

    def generate_transformation(decision_vector, centroid):
        tran = decision_vector[:3]
        rot = pso_registration.spherical_to_rotvec(*decision_vector[3:])

        rot_matrix = numpy.eye(4, 4)
        rot_matrix[:3, :3] = rot.as_dcm()
        centroid_trans = numpy.eye(4, 4)
        centroid_trans[:3, 3] = -centroid
        trans_matrix = numpy.eye(4, 4)
        trans_matrix[:3, 3] = tran
        matrix = (
            numpy.linalg.inv(centroid_trans)
            @ trans_matrix
            @ rot_matrix
            @ centroid_trans
        )
        return matrix

    def fitness(self, x):
        trans = pso_registration.generate_transformation(x, self._source_centroid)
        moved_source = copy.deepcopy(self._source)
        moved_source.transform(trans)
        # chosen = numpy.random.choice(
        #     len(moved_source.points), round(len(moved_source.points)*1)
        # )

        dist, _ = self._target_tree.query(
            numpy.asarray(moved_source.points), 1, eps=3.6, p=2, n_jobs=-1
        )
        # weights = scipy.stats.norm.pdf(dist, loc = 0, scale = 0.5)
        # weights = numpy.linalg.norm(numpy.asarray(self._target.points)[idx]-self._source_centroid,ord=2, axis=1)
        # print(weights)
        # dist = dist/weights
        threshold = numpy.quantile(dist, 0.5)
        filtered = dist[numpy.logical_and(dist > threshold / 3, dist < threshold * 3)]
        if len(filtered) < len(moved_source.points) * 0.1:
            return numpy.inf
        else:
            return [numpy.sum(filtered) / len(filtered)]

    def get_bounds(self):
        return self._bounds


def calculate_error(
    cloud1: open3d.geometry.PointCloud, cloud2: open3d.geometry.PointCloud
) -> float:
    assert len(cloud1.points) == len(
        cloud2.points
    ), "len(cloud1.points) != len(cloud2.points)"

    centroid, _ = cloud1.compute_mean_and_covariance()
    weights = numpy.linalg.norm(numpy.asarray(cloud1.points) - centroid, 2, axis=1)
    distances = numpy.linalg.norm(
        numpy.asarray(cloud1.points) - numpy.asarray(cloud2.points), 2, axis=1
    ) / len(weights)
    return numpy.sum(distances / weights)


def align(source, target, initial_trans, args):
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    source = open3d.io.read_point_cloud(source)
    target = open3d.io.read_point_cloud(target).voxel_down_sample(args.target_leaf)
    source_centroid, _ = source.compute_mean_and_covariance()
    target_ceontroid, _ = target.compute_mean_and_covariance()

    # [tx, ty, tz, inclination, radius, azimuth]
    bounds = (args.bounds[:6], args.bounds[6:])
    logging.info(f"Using {initial_trans} as initial transformation")
    logging.info(f"Using {bounds} as bounds")

    moved_source = copy.deepcopy(source)
    moved_source = moved_source.transform(initial_trans)

    initial_error = calculate_error(source, moved_source)
    moved_source = moved_source.voxel_down_sample(args.source_leaf).translate(
        target_ceontroid, relative=False
    )

    if args.display:
        vis = open3d.visualization.Visualizer()
        vis.create_window()
    prob = pygmo.problem(pso_registration(moved_source, target, bounds))


    pso = pygmo.pso(
        gen=1,
        variant=3,
        omega=0.7298,
        eta1=2.1,
        eta2=1.9,
        max_vel=0.25,
        neighb_type=2,
        neighb_param=2,
        memory=True,
    )
    # pso = pygmo.pso(
    #     gen=1,
    #     variant=5,
    #     omega=0.5,
    #     eta1=2.8,
    #     eta2=1.3,
    #     max_vel=1,
    #     neighb_type=1,
    #     memory=True,
    # )
    bee = pygmo.bee_colony(gen=1, limit=20)
    sade = pygmo.sade(gen=1, memory =True)
    algo = pygmo.algorithm(bee)

    # islands =[]
    # for i in range(args.swarm):
    #     islands.append(pygmo.island(algo = algo, prob = prob, size=args.part, udi=dummy_isl())

    archi = pygmo.archipelago(
        n=args.swarm,
        t=pygmo.ring(),
        udi=dummy_isl(),
        algo=algo,
        prob=prob,
        pop_size=args.part,
        r_pol=pygmo.fair_replace(rate=1),
        s_pol=pygmo.select_best(rate=1),
    )



    if args.display:
        source.paint_uniform_color([0, 1, 0])
        target.paint_uniform_color([1, 0, 0])
        moved_source.paint_uniform_color([0, 0, 1])

        vis.add_geometry(moved_source)
        vis.add_geometry(source)
        vis.add_geometry(target)
        previous_trans = numpy.eye(4, 4)

    for gen in range(args.gen):
        archi.evolve()
        archi.wait_check()

        

        logging.info(archi)
        best_idx = numpy.argmin(archi.get_champions_f())
        best = archi.get_champions_x()[best_idx]
        best_f = archi.get_champions_f()[best_idx]
        logging.info(f"gen: {gen} - {best} --> {best_f}")

        trans = pso_registration.generate_transformation(best, source_centroid)

        if args.display:
            moved_source.transform(trans @ numpy.linalg.inv(previous_trans))
            previous_trans = trans
            vis.update_geometry(moved_source)
            vis.poll_events()
            vis.update_renderer()

    best_idx = numpy.argmin(archi.get_champions_f())
    best = archi.get_champions_x()[best_idx]
    best_f = archi.get_champions_f()[best_idx]
    trans = pso_registration.generate_transformation(best, source_centroid)

    final_result = copy.deepcopy(source)
    final_result = final_result.transform(initial_trans)
    final_result = final_result.translate(target_ceontroid, relative=False)
    final_result = final_result.transform(trans)
    final_error = calculate_error(source, final_result)
    return (initial_error, final_error, trans)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSO registration")
    parser.add_argument("problem_file")
    parser.add_argument("folder")
    parser.add_argument("-sl", "--source-leaf", type=float, default=0)
    parser.add_argument("-tl", "--target-leaf", type=float, default=0)
    # parser.add_argument(
    #     "-m",
    #     "--matrix",
    #     type=float,
    #     nargs=12,
    #     default=numpy.eye(3, 4).flatten().tolist(),
    # )
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-d", "--display", action="store_true", default=False)
    parser.add_argument("-g", "--gen", type=int, default=300)
    parser.add_argument("-b", "--bounds", nargs=12, type=float)
    parser.add_argument("-p", "--part", type=int, default=50)
    parser.add_argument("-s", "--swarm", type=int, default=5)
    args = parser.parse_args()

    _, result_file = os.path.split(args.problem_file)
    result_file = result_file.replace(".txt", "_result.txt")

    with open(f"{result_file}", mode="w") as out_file:
        result_writer = csv.writer(
            out_file, delimiter=";", quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        result_writer.writerow(["#", args])
        result_writer.writerow(["id", "initial_error", "final_error", "transformation"])
        with open(sys.argv[1]) as csvfile:
            file_reader = csv.DictReader(csvfile, delimiter=" ")
            for row in tqdm.tqdm(file_reader):
                source = f"{args.folder}/{row['source']}"
                target = f"{args.folder}/{row['target']}"
                transformation = []
                for i in range(1, 13):
                    transformation.append(row[f"t{i}"])
                initial_trans = numpy.eye(4, 4)
                initial_trans[:3, :] = numpy.asarray(transformation).reshape(3, 4)

                (initial_error, final_error, trans) = align(
                    source, target, initial_trans, args
                )
                result = [
                    row["id"],
                    initial_error,
                    final_error,
                    " ".join([str(x) for x in trans.flatten().tolist()]),
                ]
                result_writer.writerow(result)
                out_file.flush()
