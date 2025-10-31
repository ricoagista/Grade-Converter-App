//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ReShade effect file
// visit facebook.com/MartyMcModding for news/updates
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Ambient Obscurance with Indirect Lighting "MXAO" 1.6 by Marty McFly
// Copyright © 2008-2016 Marty McFly
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

uniform float fMXAOAmbientOcclusionAmount <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 3.00;
	ui_tooltip = "MXAO: Linearly increases AO intensity. Can cause pitch black clipping if set too high.";
> = 2.00;

uniform bool bMXAOIndirectLightingEnable <
	ui_tooltip = "MXAO: Enables Indirect Lighting calculation. Will cause a major fps hit.";
> = false;

uniform float fMXAOIndirectLightingAmount <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 12.00;
	ui_tooltip = "MXAO: Linearly increases IL intensity. Can cause overexposured white spots if set too high.";
> = 4.00;

uniform float fMXAOIndirectLightingSaturation <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 3.00;
	ui_tooltip = "MXAO: Boosts IL saturation for more pronounced effect.";
> = 1.00;

uniform float fMXAOSampleRadius <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 20.00;
	ui_tooltip = "MXAO: Sample radius of GI, higher values drop performance.\nHeavily depending on game, GTASA: 2 = GTA V: 10ish.";
> = 2.50;

uniform int iMXAOSampleCount <
	ui_type = "drag";
	ui_min = 12; ui_max = 255;
	ui_tooltip = "MXAO: Amount of MXAO samples. Higher means more accurate and less noisy AO at the cost of fps.";
> = 32;

uniform bool bMXAOSmartSamplingEnable <
	ui_tooltip = "MXAO: Enables smart sample count reduction for far areas.\nEffect is lowered for low sample counts to prevent single digit sample counts in far areas.";
> = true;

uniform int iMXAOBayerDitherLevel <
	ui_type = "drag";
	ui_min = 2; ui_max = 9;
	ui_tooltip = "MXAO: 2^ditherlevel: size of AO sampling pattern size. Lower values mean less distinctive sample dirs and haloing.";
> = 7;

uniform float fMXAONormalBias <
	ui_type = "drag";
	ui_min = 0.0; ui_max = 0.8;
	ui_tooltip = "MXAO: Normals bias to reduce self-occlusion of surfaces that have a low angle to each other.";
> = 0.8;

uniform bool bMXAOPerPixelNormalsEnable <
	ui_tooltip = "MXAO: TEST! Enables per pixel normals derived from color input so surfaces get some relief instead of being 100% flat.";
> = true;

uniform bool bMXAOBackfaceCheckEnable <
	ui_tooltip = "MXAO: For indirect lighting only!\nEnables back face check so surfaces facing away from the source position don't cast light. \nIt comes with a slight fps drop.";
> = true;

uniform float fMXAOBlurSharpness <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 5.00;
	ui_tooltip = "MXAO: AO sharpness, higher means sharper geometry edges but noisier AO, less means smoother AO but blurry in the distance.";
> = 1.00;

uniform int fMXAOBlurSteps <
	ui_type = "drag";
	ui_min = 0; ui_max = 5;
	ui_tooltip = "MXAO: Offset count for AO bilateral blur filter. Higher means smoother but also blurrier AO.";
> = 3;

uniform bool bMXAODebugViewEnable <
	ui_tooltip = "MXAO: Enables raw AO/IL output for debugging and tuning purposes.";
> = false;

uniform float fMXAOFadeStart <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 1.00;
	ui_tooltip = "MXAO: Depth at which AO starts to fade out. 0.0 = camera, 1.0 = sky. Must be lower than AO fade end.";
> = 0.20;

uniform float fMXAOFadeEnd <
	ui_type = "drag";
	ui_min = 0.00; ui_max = 1.00;
	ui_tooltip = "MXAO: Depth at which AO completely fades out. 0.0 = camera, 1.0 = sky. Must be higher than AO fade start.";
> = 0.40;

#define fMXAOMipLevelIL				2	//[0 to 4]       Miplevel of IL texture. 0 = fullscreen, 1 = 1/2 screen width/height, 2 = 1/4 screen width/height and so forth. 
#define fMXAOMipLevelAO				0	//[0 to 2]	 Miplevel of AO texture. 0 = fullscreen, 1 = 1/2 screen width/height, 2 = 1/4 screen width/height and so forth. Best results: IL MipLevel = AO MipLevel + 2
#define bMXAOBoundaryCheckEnable		1	//[0 or 1]	 Enables screen boundary check for samples. Can be useful to remove odd behaviour with too high sample radius / objects very close to camera. It comes with a slight fps drop.

//custom variables, depleted after Framework implementation.
#define AO_FADE____START 		fMXAOFadeStart		//[0.0 to 1.0]	 Depth at which AO starts to fade out. 0.0 = camera, 1.0 = sky. Must be lower than AO fade end.
#define AO_FADE____END   		fMXAOFadeEnd		//[0.0 to 1.0]	 Depth at which AO completely fades out. 0.0 = camera, 1.0 = sky. Must be higher than AO fade start.

#include "ReShade.fxh"


//textures
texture2D texLOD 	{ Width = BUFFER_WIDTH; 			  Height = BUFFER_HEIGHT; 			    Format = RGBA8; MipLevels = 5+fMXAOMipLevelIL;};
texture2D texDepthLOD 	{ Width = BUFFER_WIDTH; 			  Height = BUFFER_HEIGHT;  			    Format = R16F;  MipLevels = 5+fMXAOMipLevelAO;};
texture2D texzzzzz 	{ Width = BUFFER_WIDTH; 			  Height = BUFFER_HEIGHT; 			    Format = RGBA8; MipLevels = 5+fMXAOMipLevelIL;};
texture2D texSSAO	{ Width = BUFFER_WIDTH; 	          	  Height = BUFFER_HEIGHT;            		    Format = RGBA8; 				  };

sampler2D SamplerLOD
{
	Texture = texLOD;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Clamp;
	AddressV = Clamp;
};

sampler2D SamplerDepthLOD
{
	Texture = texDepthLOD;
	MinFilter = POINT;
	MagFilter = POINT;
	MipFilter = POINT;
	AddressU = Clamp;
	AddressV = Clamp;
};

sampler2D Samplerzzzzz
{
	Texture = texzzzzz;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Clamp;
	AddressV = Clamp;
};

sampler2D SamplerSSAO
{
	Texture = texSSAO;
	MinFilter = LINEAR;
	MagFilter = LINEAR;
	MipFilter = LINEAR;
	AddressU = Clamp;
	AddressV = Clamp;
};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

float GetLinearDepth(float2 coords)
{
	return ReShade::GetLinearizedDepth(coords);
}

float3 GetPosition(float2 coords)
{
	float EyeDepth = GetLinearDepth(coords.xy)*RESHADE_DEPTH_LINEARIZATION_FAR_PLANE;
	return float3((coords.xy * 2.0 - 1.0)*EyeDepth,EyeDepth);
}

float3 GetPositionLOD(float2 coords, int mipLevel)
{
	float EyeDepth = tex2Dlod(SamplerDepthLOD, float4(coords.xy,0,mipLevel)).x;
	return float3((coords.xy * 2.0 - 1.0)*EyeDepth,EyeDepth);
}

float3 GetNormalFromDepth(float2 coords) 
{
	float3 centerPos = GetPosition(coords.xy);
	float2 offs = ReShade::PixelSize.xy;
	float3 ddx1 = GetPosition(coords.xy + float2(offs.x, 0)) - centerPos;
	float3 ddx2 = centerPos - GetPosition(coords.xy + float2(-offs.x, 0));

	float3 ddy1 = GetPosition(coords.xy + float2(0, offs.y)) - centerPos;
	float3 ddy2 = centerPos - GetPosition(coords.xy + float2(0, -offs.y));

	ddx1 = lerp(ddx1, ddx2, abs(ddx1.z) > abs(ddx2.z));
	ddy1 = lerp(ddy1, ddy2, abs(ddy1.z) > abs(ddy2.z));

	float3 normal = cross(ddy1, ddx1);
	
	return normalize(normal);
}

float3 GetNormalFromColor(float2 coords, float3 normal, float2 offset, float scale, float sharpness)
{
	float fresnelFactor = saturate(dot(normalize(GetPosition(coords.xy)),-normal.xyz));

	float2 dir1 = normalize(normal.xy)*fresnelFactor;
	float2 dir2 = normalize(float2(-normal.y, normal.x));

	float hpx = dot(tex2Dlod(ReShade::BackBuffer, float4(coords.xy + dir1 * offset,0,0)).xyz,float3(0.299,0.587,0.114)) * scale;
    	float hmx = dot(tex2Dlod(ReShade::BackBuffer, float4(coords.xy - dir1 * offset,0,0)).xyz,float3(0.299,0.587,0.114)) * scale;
    	float hpy = dot(tex2Dlod(ReShade::BackBuffer, float4(coords.xy + dir2 * offset,0,0)).xyz,float3(0.299,0.587,0.114)) * scale;
    	float hmy = dot(tex2Dlod(ReShade::BackBuffer, float4(coords.xy - dir2 * offset,0,0)).xyz,float3(0.299,0.587,0.114)) * scale;

   	float dpx = GetLinearDepth(coords.xy + dir1 * offset);
    	float dmx = GetLinearDepth(coords.xy - dir1 * offset);
    	float dpy = GetLinearDepth(coords.xy + dir2 * offset);
    	float dmy = GetLinearDepth(coords.xy - dir2 * offset);

	float2 xymult = float2(abs(dmx - dpx), abs(dmy - dpy)) * sharpness * fresnelFactor; 
	xymult = max(0.0, 1.0 - xymult);

	float ddx = (hmx - hpx) / (2.0 * length(dir1 * offset)) * xymult.x;
    	float ddy = (hmy - hpy) / (2.0 * length(dir2 * offset)) * xymult.y;
    
    	return normalize(float3(ddx, ddy, 1.0));
}

float3 GetBlendedNormals(float3 n1, float3 n2)
{
	n1 += float3( 0, 0, 1); 
	n2 *= float3(-1, -1, 1); 
	return n1*dot(n1, n2)/n1.z - n2;
}

float4 GetBlurFactors(float2 coords)
{
	return float4(tex2Dlod(Samplerzzzzz, float4(coords.xy,0,0)).xyz*2.0-1.0,GetLinearDepth(coords.xy));
}

float GetBlurWeight(float O, float4 z, float4 z0)
{
	//GFSDK with my normal angle consideration
	const float BlurSigma 		= fMXAOBlurSteps+1.0;
	const float BlurFalloff 	= 1.0 / (2.0*BlurSigma*BlurSigma);
	      float BlurFresnelFactor  	= saturate(min(-z0.z,-z.z)); 

	//we don't need any sigma etc for normals, angle is important
	float DeltaN = (1.0 + dot(-z.xyz,z0.xyz))*10.0*fMXAOBlurSharpness;
	float DeltaZ = (z.w-z0.w)*1000.0*BlurFresnelFactor*fMXAOBlurSharpness;
	return exp2(-O*O*BlurFalloff - DeltaZ*DeltaZ - DeltaN*DeltaN);
}

float GetBayerFromCoordLevel(float2 pixelpos, int maxLevel) 
{ 
	float finalBayer = 0.0; 

	for(float i = 1-maxLevel; i<= 0; i++) 
	{ 
		float bayerSize = exp2(i); 
	        float2 bayerCoord = floor(pixelpos * bayerSize) % 2.0; 
		float bayer = 2.0 * bayerCoord.x - 4.0 * bayerCoord.x * bayerCoord.y + 3.0 * bayerCoord.y; 
		finalBayer += exp2(2.0*(i+maxLevel))* bayer; 
	}

	float finalDivisor = exp2(2.0 * maxLevel + 2.0)- 4.0; 
	return finalBayer/ finalDivisor; 
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PS_AO_Pre(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 color : SV_Target0, out float4 depth : SV_Target1, out float4 zzzzz : SV_Target2)
{
	color = tex2D(ReShade::BackBuffer, texcoord.xy);
	depth = GetLinearDepth(texcoord.xy)*RESHADE_DEPTH_LINEARIZATION_FAR_PLANE;
	zzzzz.xyz = GetNormalFromDepth(texcoord.xy).xyz;

	if(bMXAOPerPixelNormalsEnable != 0)
	{
		float3 ppnormal = GetNormalFromColor(texcoord.xy,zzzzz.xyz, 70.0 * ReShade::PixelSize.xy / depth.x, 0.005, 10000.0);
		zzzzz.xyz = GetBlendedNormals(zzzzz.xyz,ppnormal.xyz);
	}

	zzzzz.xyz = zzzzz.xyz * 0.5 + 0.5;
	zzzzz.w = GetBayerFromCoordLevel(vpos.xy,iMXAOBayerDitherLevel);
}

void PS_AO_Gen(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0)
{
	float4 normalSample = tex2D(Samplerzzzzz,texcoord.xy);
	float3 ScreenSpaceNormals = normalSample.xyz * 2.0 - 1.0;
	float3 ScreenSpacePosition = GetPositionLOD(texcoord.xy, 0);

	float scenedepth = ScreenSpacePosition.z / RESHADE_DEPTH_LINEARIZATION_FAR_PLANE;
	ScreenSpacePosition += ScreenSpaceNormals * scenedepth;

	float numSamples = iMXAOSampleCount;
	float SampleRadiusScaled  = fMXAOSampleRadius / (numSamples * ScreenSpacePosition.z * 0.6);

	float rotAngle = normalSample.w;; 
	float mipFactor = SampleRadiusScaled.x*numSamples*19.0;

	float4 AOandGI = 0.0;
	float2 currentVector = 0.0;
	sincos(6.283*rotAngle, currentVector.y, currentVector.x);
	currentVector *= SampleRadiusScaled;

	float fNegInvR2 = -1.0/(fMXAOSampleRadius*fMXAOSampleRadius);
	float Aspect = ReShade::PixelSize.y/ReShade::PixelSize.x;


	[loop]
		for (float i=1.0; i <= numSamples; i++) 
	{
		currentVector = mul(currentVector.xy,  float2x2(0.575,0.81815,-0.81815,0.575));	
		float2 currentOffset = currentVector.xy * i * float2(1.0, Aspect) + texcoord.xy;
#if(bMXAOBoundaryCheckEnable != 0)
		[branch]
		if(all(saturate(-currentOffset * currentOffset + currentOffset)))
		{
#endif
			float mipLevel = clamp((int)floor(log2(mipFactor*i)) - 3, fMXAOMipLevelAO, 5); //AO must not go beyond 5

			float3 occlVec 		= -ScreenSpacePosition + GetPositionLOD(currentOffset.xy, mipLevel);
			float occlDistanceRcp 	= rsqrt(dot(occlVec,occlVec));
 			float SurfaceAngle 	= dot(occlVec, ScreenSpaceNormals)*occlDistanceRcp; 
			
			float fAO = saturate(1.0 + fNegInvR2/occlDistanceRcp)  * saturate(SurfaceAngle - fMXAONormalBias);

			if(bMXAOIndirectLightingEnable)
			{
				float3 fIL = tex2Dlod(SamplerLOD, float4(currentOffset,0,mipLevel + fMXAOMipLevelIL)).xyz;
				if(bMXAOBackfaceCheckEnable)
				{
					float3 offsetNormals = tex2Dlod(Samplerzzzzz, float4(currentOffset,0,mipLevel + fMXAOMipLevelIL)).xyz * 2.0 - 1.0; 
					//float facingtoSource = dot(-normalize(occlVec),offsetNormals);
					float facingtoSource = dot(-occlVec,offsetNormals)*occlDistanceRcp;
					facingtoSource = smoothstep(-0.5,0.0,facingtoSource); 
					fIL *= facingtoSource;
				}
				AOandGI.w += fAO*saturate(1-dot(fIL,float3(0.299,0.587,0.114)));
				AOandGI.xyz += fIL*fAO;
			}
			else
			{
				AOandGI.w += fAO;
			}
#if(bMXAOBoundaryCheckEnable != 0)
		}
#endif
	}

	//AOandGI *= 20.0 / ((1.0-fMXAONormalBias)*numSamples*fMXAOSampleRadius); 
	AOandGI /= (0.05*(1.0-fMXAONormalBias)*numSamples*fMXAOSampleRadius); 
	res = lerp(AOandGI,float4(0.0.xxx,0.0), AO_FADE____END < scenedepth); //AO FADEOUT
}

void PS_AO_Blur1(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0)
{
	//discrete sampling, I could either make this one *1.0 and in the loop every sample * 2.0 
	//or do it like this and save that multiplication inside the loop.
	float4 total_ao 		= tex2D(SamplerSSAO, texcoord.xy) * 0.5; 
	float total_weight 		= 0.5;
	float4 center_factor 		= GetBlurFactors(texcoord.xy);

	[loop]
	for(float x = 1.0; x <= fMXAOBlurSteps; x++)
	{
		float2 blurdir = float2(2.0 * x - 0.5,0.0);
		float2 blurcoord =  blurdir * ReShade::PixelSize.xy + texcoord.xy;

		float4 temp_ao = tex2Dlod(SamplerSSAO, float4(blurcoord,0,0));	
		float4 temp_factor = GetBlurFactors(blurcoord);

		float temp_weight = GetBlurWeight(blurdir.x, temp_factor, center_factor);
		total_ao += temp_ao * temp_weight;
		total_weight += temp_weight;

		blurcoord =  -blurdir * ReShade::PixelSize.xy + texcoord.xy;
		temp_ao = tex2Dlod(SamplerSSAO, float4(blurcoord,0,0));	
		temp_factor = GetBlurFactors(blurcoord);

		temp_weight = GetBlurWeight(blurdir.x, temp_factor, center_factor);
		total_ao += temp_ao * temp_weight;
		total_weight += temp_weight;
	}

	total_ao /= max(total_weight,1.0);
	res = total_ao;
}

void PS_AO_Blur2(float4 vpos : SV_Position, float2 texcoord : TEXCOORD, out float4 res : SV_Target0)
{
	float4 total_ao 		= tex2D(ReShade::BackBuffer, texcoord.xy) * 0.5;
	float total_weight 		= 0.5;
	float4 center_factor 		= GetBlurFactors(texcoord.xy);

	[loop]
	for(float y = 1.0; y <= fMXAOBlurSteps; y++)
	{
		float2 blurdir = float2(0.0, 2.0 * y - 0.5);
		float2 blurcoord = blurdir * ReShade::PixelSize.xy + texcoord.xy;

		float4 temp_ao = tex2Dlod(ReShade::BackBuffer, float4(blurcoord,0,0));	
		float4 temp_factor = GetBlurFactors(blurcoord);

		float temp_weight = GetBlurWeight(blurdir.y, temp_factor, center_factor);
		total_ao += temp_ao * temp_weight;
		total_weight += temp_weight;

		blurcoord = -blurdir * ReShade::PixelSize.xy + texcoord.xy;
		temp_ao = tex2Dlod(ReShade::BackBuffer, float4(blurcoord,0,0));	
		temp_factor = GetBlurFactors(blurcoord);

		temp_weight = GetBlurWeight(blurdir.y, temp_factor, center_factor);
		total_ao += temp_ao * temp_weight;
		total_weight += temp_weight;
	}

	total_ao /= max(total_weight,1.0);
	float4 mxao = saturate(total_ao);

	float scenedepth = GetLinearDepth(texcoord.xy); //might change center_factor so better fetch depth directly here.
	float4 color = max(0.0,tex2Dlod(SamplerLOD, float4(texcoord.xy,0,0))); //get miplevel 0 which was shifted prior in MipLODBias.
	float colorgray = dot(color.xyz,float3(0.299,0.587,0.114));

	mxao.xyz  = lerp(dot(mxao.xyz,float3(0.299,0.587,0.114)),mxao.xyz,fMXAOIndirectLightingSaturation) * fMXAOIndirectLightingAmount;
	mxao.w    = 1.0-pow(1.0-mxao.w, fMXAOAmbientOcclusionAmount * 2.0);

	if (!bMXAODebugViewEnable) mxao    = lerp(mxao, 0.0, pow(colorgray,2.0));//lerp(mxao, 0.0, saturate(colorgray * 2.0));

	mxao.w    = lerp(mxao.w, 0.0,smoothstep(AO_FADE____START, AO_FADE____END, scenedepth)); 			//AO FADEOUT
	mxao.xyz  = lerp(mxao.xyz,0.0,smoothstep(AO_FADE____START*0.5, AO_FADE____END*0.5, scenedepth)); 		//AO FADEOUT //IL can look really bad on far objects.

	float3 GI = mxao.w - mxao.xyz*2;
	GI = max(0.0,1-GI);

	color.xyz *= GI;
	color.xyz += mxao.xyz * 0.05; //pitch black surfaces.

	if (bMXAODebugViewEnable)
	{
		if (bMXAOIndirectLightingEnable)
		{
			color.xyz = GI*0.5;
		}
		else
		{
 			color.xyz = GI;
		}
	}	

	res = color;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

technique MXAO16
{
	pass P0
	{
		VertexShader = PostProcessVS;
		PixelShader  = PS_AO_Pre;
		RenderTarget0 = texLOD;
		RenderTarget1 = texDepthLOD;
		RenderTarget2 = texzzzzz;
	}
	pass P1
	{
		VertexShader = PostProcessVS;
		PixelShader  = PS_AO_Gen;
		RenderTarget = texSSAO;
	}
	pass P2_0
	{
		VertexShader = PostProcessVS;
		PixelShader  = PS_AO_Blur1;
	}
	pass P2
	{
		VertexShader = PostProcessVS;
		PixelShader  = PS_AO_Blur2;
	}
}