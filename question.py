import numpy as np

import openmc

cdte = openmc.Material(name="cdte")
cdte.add_nuclide('Cd112',1.)
cdte.add_nuclide('Te128',1.)

water = openmc.Material(name="h2o")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)

aluminium = openmc.Material(name="al")
aluminium.add_element('aluminium', 1.0)

cobalt = openmc.Material(name="co")
cobalt.add_nuclide('Co59',1.)

tungsten = openmc.Material(name="w")
tungsten.add_element('tungsten', 1.0)

cdte = openmc.Material(name="cdte")
cdte.add_nuclide('Cd112',1.)
cdte.add_nuclide('Te128',1.)

ceil = openmc.ZPlane(2500.)
floor = openmc.ZPlane(-2500.)

class Detector:
    def __init__(self):
        self.x = .35
        self.z = .35
        self.y = .175
        self.origin = (0,0,0)


detector = Detector()

class HyperParameters:
    def __init__(self):
        self.interelement = .4
        self.n_elements = 8.
        self.aperture = .15

hypp = HyperParameters()

class PicSize:
    def __init__(self):
        self.w = 80
        self.h = 60
        self.pixel_size = .1
        self.color_by = 'cell'
        #self.color_by = 'material'

ps = PicSize()

def vertical_section():
    p = openmc.Plot()
    p.filename = 'vertical'
    p.basis = 'xz'
    p.width = (ps.w,ps.h)
    p.origin = (0,0,0)
    p.pixels = (int(p.width[0]/ps.pixel_size), int(p.width[1]/ps.pixel_size))
    p.color_by = ps.color_by

    return p

def horizontal_section():
    p = openmc.Plot()
    p.filename = 'horizontal'
    p.width = (ps.w,ps.h)
    p.origin = (0,0,0)
    p.pixels = (int(p.width[0]/ps.pixel_size), int(p.width[1]/ps.pixel_size))
    p.color_by = ps.color_by

    return p


def detectors_tallies(cell_):
    detector_cells = []

    universes = cell_.fill.universes

    for u in universes:

        print(u[0].bounding_box)

        for k,v in u[0].cells.items():

            if v.name=="detector":
                print(v.region.bounding_box)
                detector_cells.append(v)

    return detector_cells


if __name__ == "__main__":

    d_prism = openmc.model.rectangular_prism(width=detector.x,height=detector.y,
        origin=(0.,0.,0.),boundary_type='transmission')

    _ = (hypp.interelement-detector.y)/4. + detector.y/2.

    u_prism = openmc.model.rectangular_prism(width=detector.x,height=hypp.interelement-detector.y/2.,
        origin=(0.,_,0.),boundary_type='transmission')

    l_prism = openmc.model.rectangular_prism(width=detector.x,height=hypp.interelement-detector.y/2.,
        origin=(0.,-_,0.),boundary_type='transmission')

    u_box = u_prism & -openmc.ZPlane(detector.z/2.) & +openmc.ZPlane(-detector.z/2.)
    l_box = l_prism & -openmc.ZPlane(detector.z/2.) & +openmc.ZPlane(-detector.z/2.)
    d_box = d_prism & -openmc.ZPlane(detector.z/2.) & +openmc.ZPlane(-detector.z/2.)

    d_cell = openmc.Cell(name="detector",fill=cdte,region=d_box)
    u_cell = openmc.Cell(fill=water,region=u_box)
    l_cell = openmc.Cell(fill=water,region=l_box)



    univ = openmc.Universe(name='detector')
    univ.add_cells([d_cell,u_cell,l_cell])

    lattice = openmc.RectLattice(name='detector')

    lattice.pitch = (detector.x, hypp.interelement)

    ll_y =-(hypp.interelement*hypp.n_elements)/2.
    ll_x = 5.
    lattice.lower_left = (ll_x,ll_y)

    lattice.universes = np.full((int(hypp.n_elements),1), univ)

    h_prism = openmc.model.rectangular_prism(width=detector.x,height=hypp.interelement*hypp.n_elements,
        origin=(ll_x+detector.x/2,0.,0.),boundary_type='transmission')

    h_box = h_prism & -openmc.ZPlane(detector.z/2.) & +openmc.ZPlane(-detector.z/2.)

    h_cell = openmc.Cell(fill=lattice,region=h_box)

    domain = openmc.model.rectangular_prism(width=1000,height=1000,
        origin=(0.,0.,0.),boundary_type='vacuum')

    domain_box = domain & -ceil & +floor & ~h_box

    domain_cell = openmc.Cell(region=domain_box,fill=water)

    model = openmc.model.Model()
    model.geometry = openmc.Geometry([h_cell,domain_cell])

    point = openmc.stats.Point((0, 0, 0))
    settings = openmc.Settings()
    src = openmc.Source(space=point,particle="photon")
    src.angle = openmc.stats.Isotropic()
    src.energy = openmc.stats.Discrete([1e6], [1.])

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.batches = 10
    settings.inactive = 2
    settings.particles = 100000
    settings.photon_transport = True

    settings.source = [src]

    model.settings = settings

    detector_cells = detectors_tallies(h_cell)

    print(detector_cells)

    tallies_list = []

    tally = openmc.Tally(name='flux')
    tally.filters = [openmc.CellFilter(detector_cells)]

    tally.scores = ['flux']

    tallies_list.append(tally)

    model.tallies = openmc.Tallies(tallies_list)

    model.export_to_xml()

    ph = horizontal_section()

    pv = vertical_section()

    plots = openmc.Plots([ph,pv])
    plots.export_to_xml()

    openmc.plot_geometry()
    openmc.run()
