#include <DD4hep/Detector.h>

void dd4hep_to_tgeo() {
  auto &detector = dd4hep::Detector::getInstance();
  detector.fromCompact("/home/andreas/cern/source/acts/acts/git2/thirdparty/OpenDataDetector/xml/OpenDataDetector.xml");
  detector.volumeManager();
  detector.apply("DD4hepVolumeManager", 0, nullptr);

  auto geometry = detector.world();
  auto& root = *geometry.placement().ptr();

  auto &manager = detector.manager();

  //root.SetVisibility(1);
  manager.GetVolume("Pixels")->SetVisibility(1);
  manager.GetVolume("PixelBarrel")->SetVisibility(0);
  manager.GetVolume("PixelEndcapP")->SetVisibility(1);
  manager.GetVolume("PixelEndcapN")->SetVisibility(1);
  //manager.GetVolume("BeamPipe")->SetVisibility(0);
  //manager.GetVolume("PST")->SetVisibility(0);
  //manager.GetVolume("Solenoid")->SetVisibility(0);

  //manager.SetVisLevel(10);
  //root.Draw("ogle");
  manager.SetVisLevel(2);
  manager.Export("odd.root");
}
