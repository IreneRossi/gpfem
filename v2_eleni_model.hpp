//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                MCiantia $
//   Date:                $Date:                    July 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_V2_ELENI_MODEL_H_INCLUDED )
#define  KRATOS_V2_ELENI_MODEL_H_INCLUDED

// System includes

// External includes
#include <iostream>
#include <fstream>

// Project includes
#include "custom_models/plasticity_models/structured_soil_model.hpp"
#include "custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.hpp"
#include "custom_models/plasticity_models/yield_surfaces/eleni_yield_surface.hpp"
#include "custom_models/elasticity_models/borja_model.hpp"
#include "custom_models/elasticity_models/tamagnini_model.hpp"



//***** the hardening law associated to this Model has ... variables
// 0. Plastic multiplier
// 1. Plastic Volumetric deformation
// 2. Plastic Deviatoric deformation
// 3. ps (mechanical)
// 4. pt (ageing)
// 5. pcSTAR = ps + (1+k) p_t
// 6. Plastic Volumetric deformation Absolut Value
// 7. NonLocal Plastic Vol Def
// 8. NonLocal Plastic Dev Def
// 9. NonLocal Plastic Vol Def ABS
// ... (the number now is then..., xD)

namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
   ///@{

   ///@name Kratos Globals
   ///@{

   ///@}
   ///@name Type Definitions
   ///@{

   ///@}
   ///@name  Enum's
   ///@{

   ///@}
   ///@name  Functions
   ///@{

   ///@}
   ///@name Kratos Classes
   ///@{

   /// Short class definition.
   /** Detail class definition.
    */
   class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) V2EleniModel : public StructuredSoilModel<TamagniniModel, EleniYieldSurface<GensNovaHardeningRule> >
   {
      public:

         ///@name Type Definitions
         ///@{

         //elasticity model
         //typedef BorjaModel                                     ElasticityModelType;
         typedef TamagniniModel                                     ElasticityModelType;
         typedef ElasticityModelType::Pointer                ElasticityModelPointer;

         //yield surface
         typedef GensNovaHardeningRule                             HardeningRuleType;
         typedef EleniYieldSurface<HardeningRuleType>    YieldSurfaceType;
         typedef YieldSurfaceType::Pointer                      YieldSurfacePointer;

         //base type
         typedef StructuredSoilModel<ElasticityModelType,YieldSurfaceType>  BaseType;

         //common types
         typedef BaseType::Pointer                         BaseTypePointer;
         typedef BaseType::SizeType                               SizeType;
         typedef BaseType::VoigtIndexType                   VoigtIndexType;
         typedef BaseType::MatrixType                           MatrixType;
         typedef BaseType::ModelDataType                     ModelDataType;
         typedef BaseType::MaterialDataType               MaterialDataType;
         typedef BaseType::PlasticDataType                 PlasticDataType;
         typedef BaseType::InternalVariablesType     InternalVariablesType;


         /// Pointer definition of V2EleniModel
         KRATOS_CLASS_POINTER_DEFINITION( V2EleniModel );

         ///@}
         ///@name Life Cycle
         ///@{

         /// Default constructor.
         V2EleniModel() : BaseType() { mInitialized = false; }

         /// Copy constructor.
         V2EleniModel(V2EleniModel const& rOther) : BaseType(rOther)  {}

         /// Assignment operator.
         V2EleniModel& operator=(V2EleniModel const& rOther)
         {
            BaseType::operator=(rOther);
            return *this;
         }

         /// Clone.
         ConstitutiveModel::Pointer Clone() const override
         {
            return ( V2EleniModel::Pointer(new V2EleniModel(*this)) );
         }

         /// Destructor.
         virtual ~V2EleniModel() {}


         ///@}
         ///@name Operators
         ///@{


         ///@}
         ///@name Operations
         ///@{


         ///@}
         ///@name Access
         ///@{


         ///@}
         ///@name Inquiry
         ///@{


         ///@}
         ///@name Input and output
         ///@{

         /// Turn back information as a string.
         virtual std::string Info() const override
         {
            std::stringstream buffer;
            buffer << "V2EleniModel" ;
            return buffer.str();
         }

         /// Print information about this object.
         virtual void PrintInfo(std::ostream& rOStream) const override
         {
            rOStream << "V2EleniModel";
         }

         /// Print object's data.
         virtual void PrintData(std::ostream& rOStream) const override
         {
            rOStream << "V2EleniModel Data";
         }


         ///@}
         ///@name Friends
         ///@{


         ///@}

      protected:
         ///@name Protected static Member Variables
         ///@{


         ///@}
         ///@name Protected member Variables
         ///@{


         ///@}
         ///@name Protected Operators
         ///@{


         ///@}
         ///@name Protected Operations
         ///@{
            // Calculate Stress and constitutive tensor
            void CalculateStressAndConstitutiveTensors(ModelDataType& rValues, MatrixType& rStressMatrix, Matrix& rConstitutiveMatrix) override
            {
               KRATOS_TRY

               // integrate "analytically" ps and pt from plastic variables. Then update the internal variables.

               PlasticDataType Variables;
               this->InitializeVariables( rValues, Variables);

               const ModelDataType & rModelData = Variables.GetModelData();
               const Properties & rMaterialProperties = rModelData.GetProperties();

               const double & rPs0 = rMaterialProperties[PS];
               const double & rPt0 = rMaterialProperties[PT];
               const double & rChis = rMaterialProperties[CHIS];
               const double & rChit = rMaterialProperties[CHIT];
               const double & rhos = rMaterialProperties[RHOS];
               const double & rhot = rMaterialProperties[RHOT];
               const double & k = rMaterialProperties[KSIM];

               const double & rPlasticVolDef = Variables.Internal.Variables[1]; 
               const double & rPlasticDevDef = Variables.Internal.Variables[2];
               const double & rPlasticVolDefAbs = Variables.Internal.Variables[6];

               double sq2_3 = sqrt(2.0/3.0);

               double ps;
               ps = rPlasticVolDef + sq2_3 * rChis * rPlasticDevDef; 
               ps = (-rPs0) * std::exp( -rhos*ps);

               double pt;
               pt = rPlasticVolDefAbs + sq2_3 * rChit * rPlasticDevDef; 
               pt = (-rPt0) * std::exp( rhot*pt);

               double pm;
               pm = ps + (1.0+k)*pt;


               mInternal.Variables[3] = ps;
               mInternal.Variables[4] = pt;
               mInternal.Variables[5] = pm;

               StructuredSoilModel::CalculateStressAndConstitutiveTensors( rValues, rStressMatrix, rConstitutiveMatrix);

               KRATOS_CATCH("")
            }
         ///@}
         ///@name Protected  Access
         ///@{


         ///@}
         ///@name Protected Inquiry
         ///@{


         ///@}
         ///@name Protected LifeCycle
         ///@{


         ///@}

      private:
         ///@name Static Member Variables
         ///@{


         ///@}
         ///@name Member Variables
         ///@{


         ///@}
         ///@name Private Operators
         ///@{


         ///@}
         ///@name Private Operations
         ///@{


         ///@}
         ///@name Private  Access
         ///@{


         ///@}
         ///@name Private Inquiry
         ///@{


         ///@}
         ///@name Serialization
         ///@{    
         friend class Serializer;

         virtual void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
         }

         virtual void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
         }

         ///@}
         ///@name Un accessible methods
         ///@{


         ///@}

   }; // Class V2EleniModel

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{


   ///@} 
   ///@name Input and output 
   ///@{


   ///@}

   ///@} addtogroup block


}  // namespace Kratos.

#endif // KRATOS_V2_ELENI_MODEL_H_INCLUDED  defined 


