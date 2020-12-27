//
// Created by Eloi on 12/19/20.
//

#ifndef DGTCM_PHYSICALMODEL_H
#define DGTCM_PHYSICALMODEL_H

namespace DGTCM {
    class PhysicalModel {
    public:
        PhysicalModel(double lengthX, double lengthY, double lengthZ,  unsigned int numXElements,  unsigned int numYElements,  unsigned int numZElements, double initialTemperature, double fluidTemperature, double rho,  double beta):
                _LengthX(lengthX), _LengthY(lengthY), _LengthZ(lengthZ),
                _NumXElements(numXElements), _NumYElements(numYElements), _NumZElements(numZElements),
                _initialTemperature(initialTemperature), _fluidTemperature(fluidTemperature),
               _rho(rho), _beta(beta) {
            _ElementXDistance = _LengthX / _NumXElements;
            _ElementYDistance = _LengthY / _NumYElements;
            _ElementZDistance = _LengthZ / _NumZElements;

            _ElementVolume = _LengthX * _LengthY * _LengthZ / (_NumXElements * _NumYElements * _NumZElements);
            _ElementXArea = _LengthY * _LengthZ / (_NumYElements * _NumZElements);
            _ElementYArea = _LengthX * _LengthZ / (_NumXElements * _NumZElements);
            _ElementZArea = _LengthX * _LengthY / (_NumXElements * _NumYElements);
            _ElementXDistance = _LengthX / _NumXElements;
            _ElementYDistance = _LengthY / _NumYElements;
            _ElementZDistance = _LengthZ / _NumZElements;
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
        double fluidTemp() const {  return _fluidTemperature; }
        double beta() const { return _beta; }
        double rho() const { return _rho; }


    private:
        double _initialTemperature;
        double _fluidTemperature;
        double _rho;
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
