#include "Camera.h"

Camera::Camera() {
	setFOV(0.25f * XM_PI, float(SCREEN_WIDTH) / (float)SCREEN_HEIGHT, 0.1f, 1000.0f);
	updateProjMat();
	updateViewMat();
}

void Camera::setFOV(float p_FOVY, float p_aspect, float p_nearPlane, float p_farPlane) {
	FOVY = p_FOVY;
	aspect = p_aspect;
	nearPlane = p_nearPlane;
	farPlane = p_farPlane;
}

void Camera::updateAspect(float p_aspect) {
	aspect = p_aspect;
	updateProjMat();
}

void Camera::rotateOnX(float angle) {
	// If rotation is between -90 and 90 degrees, rotate
	if (rotateX + angle < XM_PIDIV2 && rotateX + angle > -XM_PIDIV2) {
		rotateX += angle;
	}
}

void Camera::rotateOnY(float angle) {
	rotateY += angle;
}

void Camera::rotate() {

	// Create a rotation matrix based on the rotation angles
	XMMATRIX rotMat = XMMatrixRotationRollPitchYaw(rotateX, rotateY, 0.0f);
	// Set Right to first row of the matrix, Up to second row of the matrix, and Forward to third row of the matrix
	XMStoreFloat3(&right, rotMat.r[0]);
	XMStoreFloat3(&up, rotMat.r[1]);
	XMStoreFloat3(&forward, rotMat.r[2]);
}

void Camera::translate(XMFLOAT3 distance) {
	XMVECTOR P = XMLoadFloat3(&position);
	XMVECTOR D = XMLoadFloat3(&distance);
	D *= MOVE_SCALAR;
	XMStoreFloat3(&position, XMVectorAdd(P, D));
}

void Camera::updateViewMat() {
	// Use XMMatrixLookToLH to create a view matrix
	XMMATRIX V = XMMatrixLookToLH(XMLoadFloat3(&position), XMLoadFloat3(&forward), XMLoadFloat3(&up));
	XMStoreFloat4x4(&viewMat, V);
	XMStoreFloat4x4(&viewProjMat, XMLoadFloat4x4(&projMat) * XMLoadFloat4x4(&viewMat));
}

void Camera::updateProjMat() {
	XMMATRIX P = XMMatrixPerspectiveFovLH(FOVY, aspect, nearPlane, farPlane);
	XMStoreFloat4x4(&projMat, P);

	XMStoreFloat4x4(&viewProjMat, XMLoadFloat4x4(&projMat) * XMLoadFloat4x4(&viewMat));
}

XMMATRIX Camera::getViewMat() {
	return XMLoadFloat4x4(&viewMat);
}

XMMATRIX Camera::getProjMat() {
	return XMLoadFloat4x4(&projMat);
}

XMMATRIX Camera::getViewProjMat() {
	return XMLoadFloat4x4(&viewProjMat);
}

XMMATRIX Camera::getInvViewProjMat()
{
	return XMMatrixInverse(nullptr, XMLoadFloat4x4(&viewProjMat));
}

void Camera::kmStateCheck(DirectX::Keyboard::State kState, DirectX::Mouse::State mState) {
	
	if (kState.W) {
		// Translate on forward
		translate({ forward.x, 0, forward.z });
	}
	if (kState.A) {
		// Translate on Left
		translate({ -right.x, 0, -right.z });
	}
	if (kState.S) {
		// Translate on Back
		translate({ -forward.x, 0, -forward.z });
	}
	if (kState.D) {
		translate({ right.x, 0, right.z });
	}
	if (kState.Space) {
		translate({ 0.f, 1.0f, 0.f });
	}
	if (kState.LeftControl) {
		translate({ 0.f, -1.0f, 0.f });
	}

	if (mState.positionMode == Mouse::MODE_RELATIVE && kState.LeftShift) {
		rotateOnX(mState.y * 0.01f);
		rotateOnY(mState.x * 0.01f);
		rotate();
	}

	updateViewMat();
}
