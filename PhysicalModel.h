//
// Created by Eloi on 12/19/20.
//

#ifndef DGTCM_PHYSICALMODEL_H
#define DGTCM_PHYSICALMODEL_H

namespace DGTCM {
    class PhysicalModel {
    public:
        PhysicalModel(double lengthX, double lengthY, double lengthZ,  unsigned int numXElements,  unsigned int numYElements,  unsigned int numZElements, double initialTemperature, double fluidTemperature, double alpha, double rho, double Cp, double beta):
                _LengthX(lengthX), _LengthY(lengthY), _LengthZ(lengthZ),
                _NumXElements(numXElements), _NumYElements(numYElements), _NumZElements(numZElements),
                _initialTemperature(initialTemperature), _fluidTemperature(fluidTemperature),
                _alpha(alpha),_rho(rho), _Cp(Cp), _beta(beta) {
            _ElementXDistance = lengthX / numXElements;
            _ElementYDistance = lengthY / numYElements;
            _ElementZDistance = lengthZ / numZElements;

            _ElementVolume = lengthX * lengthY * lengthZ / (numXElements * numYElements * numZElements);
            _ElementXArea = lengthY * lengthZ / (numYElements * numZElements);
            _ElementYArea = lengthX * lengthZ / (numXElements * numZElements);
            _ElementZArea = lengthX * lengthY / (numXElements * numYElements);
            _ElementXDistance = lengthX / numXElements;
            _ElementYDistance = lengthY / numYElements;
            _ElementZDistance = lengthZ / numZElements;
        }

        double nX() const { return _NumXElements; }
        double nY() const { return _NumYElements; }
        double nZ() const { return _NumZElements; }
        double vol() const { return _ElementVolume; }
        double xArea() const { return _ElementXArea; }
        double yArea() const { return _ElementYArea; }
        double zArea() const { return _ElementZArea; }
        double dX() const { return _ElementXDistance; }
        double dY() const { return _ElementYDistance; }
        double dZ() const { return _ElementZDistance; }
        double initialTemp() const { return _initialTemperature; }
        double alpha() const { return _alpha; }
        double fluidTemp() const {  return _fluidTemperature; }
        double beta() const { return _beta; }
        double rho() const { return _rho; }
        double Cp() const { return _Cp; }


    private:
        double _initialTemperature;
        double _fluidTemperature;
        double _alpha;
        double _rho;
        double _Cp;
        double _beta;
        double _LengthX;
        double _LengthY;
        double _LengthZ;
        unsigned int _NumXElements;
        unsigned int _NumYElements;
        unsigned int _NumZElements;
        double _ElementVolume;
        double _ElementXArea;
        double _ElementYArea;
        double _ElementZArea;
        double _ElementXDistance;
        double _ElementYDistance;
        double _ElementZDistance;
    };
}


#endif //DGTCM_PHYSICALMODEL_H
