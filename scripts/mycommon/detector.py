import acts
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector


u = acts.UnitConstants

# ODD configs
geoDir = getOpenDataDetectorDirectory()
materialMap = geoDir / "data/odd-material-maps.root"
digiConfig = geoDir / "config/odd-digi-smearing-config.json"
seedingSel = geoDir / "config/odd-seeding-config.json"
materialDeco = acts.IMaterialDecorator.fromFile(materialMap)

# ODD
detector, trackingGeometry, decorators = getOpenDataDetector(
    geoDir, mdecorator=materialDeco
)
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))


def get_odd():
    return detector, trackingGeometry, decorators, field, digiConfig, seedingSel
