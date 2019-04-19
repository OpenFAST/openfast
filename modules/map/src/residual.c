/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   map[dot]plus[dot]plus[dot]help[at]gmail                    *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include "residual.h"


double residual_function_length_no_contact(const double V, const double H, const double w, const double Lu, const double EA, const double l)
{  
  return (H/w)*ARCSINH(V/H) - (H/w)*ARCSINH( (V-w*Lu)/H ) + ((H*Lu)/(EA)) - l;
};


double residual_function_height_no_contact(const double V, const double H, const double w, const double Lu, const double EA, const double h)
{
  return (H/w)* sqrt(1 + pow((V/H), 2)) - (H/w)*sqrt(1 + pow(((V-w*Lu)/H), 2)) + 1/(EA)*(V*Lu - (w*Lu*Lu)/2) - h;
};


double residual_function_length_contact(const double V, const double H, const double w, const double Lu, const double EA, const double l, const double cb)
{  
  /* Note that Lb = Lu - V/w */
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*(Lu-V/w)*(Lu-V/w) + (Lu/EA)*H + (Lu-V/w) - l;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return log((V/H) + sqrt(1.0 + pow(V/H,2)))*(H/w) - 0.5*(cb/EA)*w*((Lu-V/w)*(Lu-V/w) - ((Lu-V/w) - (H/w)/cb)*((Lu-V/w) - (H/w)/cb)) + (Lu/EA)*H + (Lu-V/w) - l;
  };
};


double residual_function_height_contact(const double V, const double H, const double w, const double Lu, const double EA, const double h, const double cb)
{
  // @todo remove this conditional statement because the equations are the same regardless of the outcome
  if (-cb*(V-w*Lu)<H) { /* true when a portion of the line rests on the seabed and the anchor tension is nonzero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  } else { /* 0.0<H<=-CB*(V-w*Lu), A  portion of the line must rest on the seabed and the anchor tension is zero */
    return (sqrt(1.0 + pow(V/H,2)) - 1.0)*(H/w) + 0.5*V*(V/(w*EA)) - h;
  };
};
