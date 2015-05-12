#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "pbrt.h"
#include "camera.h"
#include "film.h"
#include "disk.h"
#include "sphere.h"
#include "diffgeom.h"

class AfZone {
	public:
	  float left, right;
	  float top, bottom;
	  int xres,yres;
};

class Lens {
public:
    float radius;
    float zPos;
    float nReflect;
    float aperture;
    
    Lens(float r, float zI, float n, float ap):
        radius(r), zPos(zI), nReflect(n), aperture(ap) {};
};

class LensMask {
public:
    float radius;
    float zPos;
    float aperture;
    float surfArea;
    
    LensMask() {};
    LensMask(Lens& lastLens) {
        radius = lastLens.radius;
        zPos = lastLens.zPos;
        aperture = lastLens.aperture;
        surfArea = pow(lastLens.aperture / 2.f, 2) * M_PI;
    };
    Point getSamplePoint(const CameraSample &sample) const;
};

class RealisticCamera : public Camera {
public:
    RealisticCamera(const AnimatedTransform &cam2world,
      float hither, float yon, float sopen,
      float sclose, float filmdistance, float aperture_diameter,
      const string &specfile,
      const string &autofocusfile,
      float filmdiag,
      Film *film);
    ~RealisticCamera();
    float GenerateRay(const CameraSample &sample, Ray *) const;
    void  AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample);
    void ParseAfZones(const string& filename);
    void ParseSpec(const string& filename);

private:
    bool  autofocus;
    vector<AfZone> afZones;
    vector<Lens> lenses;
    float ShutterOpen;
    float ShutterClose;
    Film * film;
    
    float filmZPos;
    float filmDiag;
    float hitherP;
    float yonP;
    Transform raster2camera;
    
    int apertureIdx;
    LensMask lensMask;
    
    bool traceThruLenses(Ray* ray) const;
    float getFocusMeasure(AfZone & zone, Renderer * renderer, const Scene * scene, Sample * origSample);
    float computeScoreFromRGB(float *rgb, float *lumin, int width, int height);
    float estimateDepth(vector<float> filmPosition, vector<float> measureScore);
    void RGBArray2Luminance(const float *rgb, float** lumin, int width, int height);
    void computeRaster2CameraTransform();
};

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);

#endif
