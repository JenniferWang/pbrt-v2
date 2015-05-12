// cameras/realistic.cpp*
#include "stdafx.h"
#include "cameras/realistic.h"
#include "paramset.h"
#include "sampler.h"
#include "montecarlo.h"
#include "filters/box.h"
#include "film/image.h"
#include "samplers/stratified.h"
#include "intersection.h"
#include "renderer.h"

#include <fstream>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

static inline int max(int i, int j) {
    return i < j ? j : i;
}

static inline int min(int i, int j) {
    return i < j ? i : j;
}

Point LensMask::getSamplePoint(const CameraSample &sample) const{
 
    float lensU, lensV, lensZ, z;
    ConcentricSampleDisk(sample.lensU,sample.lensV,&lensU,&lensV);
    lensU *= aperture / 2.f;
    lensV *= aperture / 2.f;
    
    z = sqrtf(radius * radius - lensV * lensV - lensU * lensU);
    lensZ = zPos - radius - ((radius < 0) ? z : -z);
    
    return Point(lensU, lensV, lensZ);
}


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	   // Extract common camera parameters from \use{ParamSet}
	   float hither = params.FindOneFloat("hither", -1);
	   float yon = params.FindOneFloat("yon", -1);
	   float shutteropen = params.FindOneFloat("shutteropen", -1);
	   float shutterclose = params.FindOneFloat("shutterclose", -1);

	   // Realistic camera-specific parameters
	   string specfile = params.FindOneString("specfile", "");
	   float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
	   float fstop = params.FindOneFloat("aperture_diameter", 1.0);
	   float filmdiag = params.FindOneFloat("filmdiag", 35.0);
	   string autofocusfile = params.FindOneString("af_zones", "");
	   assert(hither != -1 && yon != -1 && shutteropen != -1 &&
	      shutterclose != -1 && filmdistance!= -1);
	   if (specfile == "") {
	       Severe( "No lens spec file supplied!\n" );
	   }
	   return new RealisticCamera(cam2world, hither, yon,
	      shutteropen, shutterclose, filmdistance, fstop,
	      specfile, autofocusfile, filmdiag, film);
}


void RealisticCamera::computeRaster2CameraTransform() {
    float diag = sqrtf( film->xResolution * film->xResolution + film->yResolution * film->yResolution );
    float ratio = filmDiag / diag;
    float xScreen = ratio * 0.5 * film->xResolution;
    float yScreen = ratio * 0.5 * film->yResolution;
    raster2camera = Translate(Vector(0.f, 0.f, filmZPos)) *
    Translate(Vector(xScreen, -yScreen, 0.f)) *
    Scale(ratio, ratio, 1) * Scale(-1.f, 1.f, 1.f);
    
}


RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
                                 float hither, float yon,
                                 float sopen, float sclose,
                                 float filmdistance, float aperture_diameter_,
                                 const string &specfile,
								 const string &autofocusfile,
                                 float filmdiag,
								 Film *f)
                                 : Camera(cam2world, sopen, sclose, f),
								   ShutterOpen(sopen),
								   ShutterClose(sclose),
								   film(f)
{

    ParseSpec(specfile);
    yonP = yon;
    filmDiag = filmdiag;
    filmZPos = lenses.back().zPos - filmdistance;
    hitherP = hither;
    
    computeRaster2CameraTransform();
    
    // update apertureStop in lenses
    lenses[apertureIdx].aperture = aperture_diameter_;
    lenses[apertureIdx].nReflect = 1.f;

    // create lensMask
    lensMask = LensMask(lenses.back());
    
	// If 'autofocusfile' is the empty string, then you should do
	// nothing in any subsequent call to AutoFocus()
	autofocus = false;

	if (autofocusfile.compare("") != 0)  {
		ParseAfZones(autofocusfile);
		autofocus = true;
	}
}


void RealisticCamera::ParseAfZones(const string& filename) {
    
  ifstream specfile(filename.c_str());
   if (!specfile) {
      fprintf(stderr, "Cannot open file %s\n", filename.c_str());
      exit (-1);
   }

   char line[512];

   while (!specfile.eof()) {
      specfile.getline(line, 512);
      if (line[0] != '\0' && line[0] != '#' &&
         line[0] != ' ' && line[0] != '\t' && line[0] != '\n') {
          
		afZones.resize(afZones.size() + 1);
		AfZone& zone = afZones[afZones.size()-1];
		sscanf(line, "%f %f %f %f\n", &zone.left, &zone.right, &zone.top, &zone.bottom);
      }
   }

	printf("Read in %zu AF zones from %s\n", afZones.size(), filename.c_str());
}


void RealisticCamera::ParseSpec(const string& filename) {
    ifstream specfile(filename.c_str());
    if (!specfile) {
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        exit (-1);
    }
    
    char line[512];
    
    float dist = 0;
    float r, zSep, nReflect, aperture;
    while (!specfile.eof()) {
        specfile.getline(line, 512);
        if (line[0] != '\0' && line[0] != '#' &&
            line[0] != ' ' && line[0] != '\t' && line[0] != '\n') {
            
            sscanf(line, "%f %f %f %f\n", &r, &zSep, &nReflect, &aperture);
            Lens lens = Lens(r, dist, nReflect, aperture);
            dist -= zSep;
            
            if (r == 0) {
                apertureIdx = lenses.size();
            }
            lenses.push_back(lens);
        }
    }
    printf("Read in %zu lenses from %s\n", lenses.size(), filename.c_str());
}


RealisticCamera::~RealisticCamera() { }


static bool traceThruOneLens(Ray& ray, const Lens &lens,
                             const float enterReflec, const float exitReflec) {
    
    // compute ray-sphere intersection
    // solve the quadratic equation ad^2 + bd + c = 0;
    // a = ||ray.d|| = 1;
    Point lensCenter(0.f, 0.f, lens.zPos - lens.radius);
    
    if (lens.radius < 0 && (ray.o - lensCenter).LengthSquared() < lens.radius * lens.radius ) {
        return false;
    }
    
    float a = 1.f;
    float b = 2 * Dot(ray.d, ray.o - lensCenter);
    float c = (ray.o - lensCenter).LengthSquared() - lens.radius * lens.radius;
    float delta = b * b - 4 * a * c;
    float t;
    
    if (delta < 0) {
        return false;
    }
    else if (delta == 0) {
        t = - b / 2.f;
    }
    else {
        if (b > 0 && lens.radius < 0)
            return false;
        
        if (lens.radius < 0)
            t = ( -b - sqrt(delta) ) / 2.f;
        else
            t = ( -b + sqrt(delta) ) / 2.f;
        
        if (t < 0)
            return false;
    }
    // update ray to the intersection point
    ray.o += ray.d * t;
    
    // check whether within the aperture
    if (lens.radius < 0 && ((ray.o - lensCenter).z > 0 || sqrt(lens.radius * lens.radius - (ray.o - lensCenter).LengthSquared()) > lens.aperture / 2.f )) {
        return  false;
    }
    
    if (lens.radius > 0 && ((ray.o - lensCenter).z < 0 || sqrt(lens.radius * lens.radius - (ray.o - lensCenter).LengthSquared()) > lens.aperture / 2.f )) {
        return false;
    }
    
    // use snell's law to compute the new direction
    Vector s1 = ray.d, n = lens.radius < 0 ? Normalize(ray.o - lensCenter): -Normalize(ray.o - lensCenter) ;
    float reflecRatio = enterReflec / exitReflec;
    Vector ns1 = Cross(n, s1);
    if (1 - (reflecRatio * reflecRatio) * Dot(ns1, ns1) < 0)
        return false;
    Vector s2 = reflecRatio * Cross(n, -ns1)
                - n * sqrt(1 - (reflecRatio * reflecRatio) * Dot(ns1, ns1));

    ray.d = Normalize(s2);
    return true;
}


bool RealisticCamera::traceThruLenses(Ray* ray) const {

    assert(lenses.size() > 1);
    int num = lenses.size();
    float enterReflec = lenses[num - 1].nReflect, exitReflec;

    for (int i = num - 1; i >= 0; i--) {
        exitReflec = (i == 0 ? 1 : lenses[i - 1].nReflect);
        if (i == apertureIdx) {
            float deltaZ = lenses[i].zPos - ray->o.z;
            float t = deltaZ / ray->d.z;
            ray->o += t * ray->d;
            if (ray->o.x * ray->o.x + ray->o.y * ray->o.y > lenses[i].aperture * lenses[i].aperture)
                return false;
        }
        else {
            if (!traceThruOneLens(*ray, lenses[i], enterReflec, exitReflec)) {
                return false;
            }
        }
        enterReflec = exitReflec;
    }
    return true;
}


float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
    

    Point pRas(sample.imageX, sample.imageY, 0);
    Point pCam;
    raster2camera(pRas, &pCam);
    Point pLens = lensMask.getSamplePoint(sample);
    
    // Set the ray
    ray->o = pCam;
    ray->d = Normalize(pLens-pCam);
    ray->mint = 0.f;
    ray->maxt = INFINITY;

    if(!traceThruLenses(ray)) {
        return 0.f;
    }
    else {
        ray->mint = hitherP;
        ray->maxt = ray->maxt = (yonP - hitherP) / ray->d.z;
        ray->time = Lerp(sample.time, shutterOpen, shutterClose);
        CameraToWorld(*ray, ray);
        ray->d = Normalize(ray->d);
        
        // compute weight
        Vector mask2Film = Normalize(pLens - pCam);
        float cosTheta = Dot(Vector(0.f, 0.f, 1.f), mask2Film);
        float weight = (lensMask.surfArea / pow(fabs(filmZPos - pLens.z), 2.f)) * pow(cosTheta, 4.f);
        if (weight > 0.3)
            printf("%f\n", weight);
        return weight;
    }
}


float RealisticCamera::getFocusMeasure(AfZone & zone, Renderer * renderer, const Scene * scene, Sample * origSample) {
    
    RNG rng;
    MemoryArena arena;
    Filter * filter = new BoxFilter(.5f,.5f);
    const float crop[] = {zone.left,zone.right,zone.top,zone.bottom};
    
    ImageFilm sensor(film->xResolution, film->yResolution, filter, crop,"foo.exr",false);
    
    int xstart,xend,ystart,yend;
    sensor.GetSampleExtent(&xstart,&xend,&ystart,&yend);
    
    StratifiedSampler sampler(xstart, xend, ystart, yend,
                              16, 16, true, ShutterOpen, ShutterClose);
    
    // Allocate space for samples and intersections
    int maxSamples = sampler.MaximumSampleCount();
    Sample *samples = origSample->Duplicate(maxSamples);
    RayDifferential *rays = new RayDifferential[maxSamples];
    Spectrum *Ls = new Spectrum[maxSamples];
    Spectrum *Ts = new Spectrum[maxSamples];
    Intersection *isects = new Intersection[maxSamples];
    
    // Get samples from _Sampler_ and update image
    int sampleCount;
    while ((sampleCount = sampler.GetMoreSamples(samples, rng)) > 0) {
        // Generate camera rays and compute radiance along rays
        for (int i = 0; i < sampleCount; ++i) {
            // Find camera ray for _sample[i]_
            float rayWeight = this->GenerateRayDifferential(samples[i], &rays[i]);
            rays[i].ScaleDifferentials(1.f / sqrtf(sampler.samplesPerPixel));
        
            // Evaluate radiance along camera ray
            if (rayWeight > 0.f)
                Ls[i] = rayWeight * renderer->Li(scene, rays[i], &samples[i], rng,
                                                 arena, &isects[i], &Ts[i]);
            else {
                Ls[i] = 0.f;
                Ts[i] = 1.f;
            }
            // Issue warning if unexpected radiance value returned
            if (Ls[i].HasNaNs()) {
                Error("Not-a-number radiance value returned "
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            else if (Ls[i].y() < -1e-5) {
                Error("Negative luminance value, %f, returned"
                      "for image sample.  Setting to black.", Ls[i].y());
                Ls[i] = Spectrum(0.f);
            }
            else if (isinf(Ls[i].y())) {
                Error("Infinite luminance value returned"
                      "for image sample.  Setting to black.");
                Ls[i] = Spectrum(0.f);
            }
            
        }
        
        // Report sample results to _Sampler_, add contributions to image
        if (sampler.ReportResults(samples, rays, Ls, isects, sampleCount)) {
            for (int i = 0; i < sampleCount; ++i) {
                sensor.AddSample(samples[i], Ls[i]);
            }
        }
        
        // Free _MemoryArena_ memory from computing image sample values
        arena.FreeAll();
    }
    
    float * rgb;
    int width;
    int height;
    sensor.WriteRGB(&rgb,&width,&height,1.f);

    float *lumin = new float[width * height];
    float score = computeScoreFromRGB(rgb, lumin, width, height);
    
    delete [] rgb;
    delete [] lumin;

    sensor.WriteImage(0.5);
    
    delete[] samples;
    delete[] rays;
    delete[] Ls;
    delete[] Ts;
    delete[] isects;
    
    return score;
}


void RealisticCamera::RGBArray2Luminance(const float *rgb, float** lumin, int width, int height) {
    
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float xyz[3];
            RGBToXYZ(&rgb[i * width * 3 + j * 3], xyz);
            (*lumin)[i * width + j] = xyz[1];
        }
    }
}


float RealisticCamera::computeScoreFromRGB(float *rgb, float *lumin, int width, int height) {
    
    assert(width > 0 && height > 0 && rgb != NULL);
    RGBArray2Luminance(rgb, &lumin, width, height);
    
    float score = 0, modifiedLaplacian, threshold = 0.1;
    int step = 2, N = 1;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            
            for (int x = max(0, i - N) ; x <= min(height - 1, i + N); x++) {
                for (int y = max(0, j - N); y <= min(width - 1, j + N); y++) {
                    
                    modifiedLaplacian = fabsf(2 * lumin[ x * width + y ]
                                              - lumin[ max(x - step, 0) * width + y ]
                                              - lumin[ min(x + step, height - 1) * width + y ])
                                        + fabsf(2 * lumin[ x * width + y ]
                                                - lumin[ x * width + min(y + step, width - 1)]
                                                - lumin[ x * width + max(y - step, 0)]);
                    score += modifiedLaplacian > threshold ? modifiedLaplacian : 0;
                }
            }
            
        }
    }
    return score;
}


void  RealisticCamera::AutoFocus(Renderer * renderer, const Scene * scene, Sample * origSample) {
    
    printf("The original film z pos is %f, ", filmZPos);
    
    if (!autofocus)
        return;
    int measureSetSize = 20;
    float minFilmZPos = filmZPos - 50;
    float maxFilmZPos = filmZPos + min(50, lenses.back().zPos - filmZPos - 10);
    float deltaDist = (maxFilmZPos - minFilmZPos) / measureSetSize;
    filmZPos = minFilmZPos;
    
    //Currently only test one AF zone.
    vector<float> filmPosition, measureScore;
    for (int i = 0; i < measureSetSize; i++) {
        computeRaster2CameraTransform();
        filmPosition.push_back(filmZPos);
        measureScore.push_back(getFocusMeasure(afZones.back(), renderer, scene, origSample));
        printf("filmZPos = %f, score = %f\n", filmZPos, measureScore.back());
        filmZPos += deltaDist;
    }
    filmZPos = estimateDepth(filmPosition, measureScore);
    printf("After af, z pos is now %f, ", filmZPos);
    computeRaster2CameraTransform();
}


float RealisticCamera::estimateDepth(vector<float> filmPosition, vector<float> measureScore) {
    
    // use the algorithm in Shree K. (1992)
    assert(measureScore.size() > 2);
    
    int k = 2;
    float fm[] = {0.f, 0.f, 0.f}, dm[] = {0.f, 0.f, 0.f};
    float deltaD = filmPosition[1] - filmPosition[0];
    
    while (k < measureScore.size()) {
        if (measureScore[k - 1] > fm[1]
            && measureScore[k - 1] > measureScore[k]
            && measureScore[k - 1] > measureScore[k - 2] ) {
            
            fm[1] = measureScore[k - 1];
            fm[0] = measureScore[k - 2];
            fm[2] = measureScore[k];
            dm[1] = filmPosition[k - 1];
        }
        k++;
    }
    dm[0] = dm[1] - deltaD;
    dm[2] = dm[1] + deltaD;
    
    float logFm [] = {log(fm[0]), log(fm[1]), log(fm[2])};
    float meanD = ((logFm[1] - logFm[2]) * (dm[1] * dm[1] - dm[0] * dm[0])
                   - (logFm[1] - logFm[0]) * (dm[1] * dm[1] - dm[2] * dm[2]))
                   / (2 * deltaD * ((logFm[1] - logFm[0]) + (logFm[1] - logFm[2])));
    
    if (!meanD) {
        return dm[1];
    }
    return meanD;
}

